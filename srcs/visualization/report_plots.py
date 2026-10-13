"""Render saved analysis results; this module never reruns identification/permutations."""
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import logging
from pathlib import Path
import pickle
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

from srcs.utils.io import save_json, _atomic
from srcs.utils.progress import emit_progress, progress_items
from srcs.visualization.plot_series import PlotSeries


@dataclass
class FigureJob:
    path: str
    kind: str
    data: object
    style: dict


def comparison_folder(value):
    """Use explicit comparison axes, never encoded display labels or hashes."""
    mode = value['name'].split(':', 1)[0]
    if mode not in ('observation_start', 'adjacent', 'all_pairs', 'windows'):
        raise ValueError(f'Unknown comparison mode: {mode}')
    if mode == 'windows':
        return 'report/figures/identification/windows'
    first, second = value['first'], value['second']
    if len(first) != 1 or len(second) != 1:
        raise ValueError('A week-pair comparison requires one week on each axis')
    return f'report/figures/identification/{mode}/{first[0]}_to_{second[0]}'


def lower_triangle_grid(ax, size):
    """Historical Phi/DP staircase boundary, including diagonal cells."""
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(False)
    for i in range(size):
        ax.hlines(i - .5, -.5, i + .5, color='black', linewidth=.5, clip_on=False)
        ax.vlines(i + .5, i - .5, size - .5, color='black', linewidth=.5, clip_on=False)
    ax.hlines(size - .5, -.5, size - .5, color='black', linewidth=.5, clip_on=False)
    ax.vlines(-.5, -.5, size - .5, color='black', linewidth=.5, clip_on=False)


def structure_metrics(ctx, results):
    """Historical DD/ND descriptions on unique observation pairs, not inference."""
    rows, seen = [], set()
    supplement = ctx.section('supplement')
    atlases = supplement.get('structure_atlases', ctx.atlases)
    names = supplement.get('structure_group_names', ['_'.join(a) for a in atlases])
    labels = dict(zip(('_'.join(a) for a in atlases), names))
    for key, value in results.items():
        atlas = key.split('/', 1)[0]
        if value.get('status') != 'complete' or atlas not in labels or atlas in supplement.get('exclude', []):
            continue
        pair = (atlas, tuple(value['first']), tuple(value['second']))
        if pair in seen:
            continue
        seen.add(pair)
        matrix = np.asarray(value['similarity_z'], dtype=float)
        matrix = matrix - matrix.min()
        total = matrix.sum()
        i, j = np.indices(matrix.shape)
        rows.append(dict(atlas=atlas, group=labels[atlas], comparison=value['name'],
                         first_weeks=';'.join(value['first']), second_weeks=';'.join(value['second']),
                         mice=';'.join(value['mice']),
                         DD=np.trace(matrix)/total if total else np.nan,
                         ND=np.sum(matrix*np.abs(i-j))/total/(len(matrix)-1) if total else np.nan))
    return pd.DataFrame(rows, columns=['atlas', 'group', 'comparison', 'first_weeks', 'second_weeks', 'mice', 'DD', 'ND'])


def edge_proportions(matrix, labels, threshold):
    """Original circularGraph matrix: signed > threshold; symmetric edge counts / 2E."""
    prefixes = [str(label).split('-')[0] for label in labels]
    groups = sorted(set(prefixes))
    group_index = {label: i for i, label in enumerate(groups)}
    counts = np.zeros((len(groups), len(groups)))
    i, j = np.triu_indices(len(labels), 1)
    selected = np.asarray(matrix)[i, j] > threshold
    for a, b in zip(i[selected], j[selected]):
        a, b = group_index[prefixes[a]], group_index[prefixes[b]]
        counts[a, b] += 1
        counts[b, a] += 1
    if selected.sum():
        counts /= 2*selected.sum()
    counts[np.triu_indices(len(groups), 1)] = np.nan
    return counts, groups


def build_jobs(ctx, records, results, reliability, frame):
    jobs = []
    paths = set()

    def add(path, kind, data, style=None):
        # Validate instead of silently rewriting user IDs or adding opaque suffixes.
        for part in path.split('/'):
            if (not part or part in ('.', '..') or part.endswith((' ', '.'))
                    or any(c in part for c in '<>:"\\|?*') or any(ord(c) < 32 for c in part)
                    or part.split('.')[0].upper() in {'CON', 'PRN', 'AUX', 'NUL',
                        *(f'COM{i}' for i in range(1, 10)), *(f'LPT{i}' for i in range(1, 10))}):
                raise ValueError(f'Invalid output path component: {part!r}')
        if path.casefold() in paths:
            raise ValueError(f'Duplicate figure path: {path}')
        paths.add(path.casefold())
        jobs.append(FigureJob(path, kind, data, dict(style or {})))

    for (atlas, week, mouse, run), matrix in records.items():
        stem = f'{week}_{mouse}_{run}_FC_Matrix_{atlas}_excluded'
        add(f'Correlation/Figures/{stem}.png', 'fc', matrix, ctx.viz_style('fc_heatmap'))
    add('report/figures/statistics/trajectories.png', 'trajectory', frame)

    complete = {key: value for key, value in results.items() if value.get('status') == 'complete'}
    conditions = {}
    for key, value in complete.items():
        atlas, condition = key.split('/', 1)
        folder = f'{comparison_folder(value)}/{atlas}'
        title = f'{atlas}\n{condition}'
        matrix_style = ctx.viz_style('id_similarity_matrix')
        matrix_style.update(title=title, xlabel='Second observation: mouse', ylabel='First observation: mouse')
        add(folder + '/similarity.png', 'similarity', (value['similarity_z'], value['mice']), matrix_style)
        add(folder + '/paired.png', 'paired_vertical', [value['other_z'], value['self_z']],
            dict(ctx.viz_style('id_paired'), title=title, xticklabels=['Other mice', 'Self'], annotate_stats=False))
        if value.get('permutation_accuracy') is not None:
            add(folder + '/permutation.png', 'hist', value['permutation_accuracy'],
                dict(ctx.viz_style('id_shuffle_hist'), arrow_x=value['accuracy'],
                     title=f"{title}\nPermutation p = {value['p_permutation']:.4g}"))
        labels = value['roi_labels']
        iu = np.triu_indices(len(labels), 1)
        for metric in ('phi', 'dp'):
            matrix = np.zeros((len(labels), len(labels)))
            vector = np.asarray(value['edge'][metric], dtype=float)
            if vector.shape != (len(iu[0]),):
                raise ValueError(f'Edge/ROI shape mismatch: {key}/{metric}')
            matrix[iu] = vector
            matrix += matrix.T
            opts = ctx.section('Edge')
            # Keep the historical threshold convention (percentile of the full matrix).
            style = dict(ctx.viz_style('edge_' + metric),
                         threshold=float(np.percentile(matrix, opts.get('threshold', 99))),
                         gap=opts.get('gap', .25), cross_only=opts.get('cross_only', False),
                         title=f'{title}\n{metric.upper()} (exploratory)', figsize=(8, 8))
            add(folder + f'/{metric}_circular.png', 'circular', (matrix, labels), style)
            if opts.get('save_mat', True):
                add(folder + f'/{metric}_matrix.png', 'edge_proportion', (matrix, labels),
                    dict(title=style['title'] + '\nSelected edge proportions',
                         threshold=style['threshold'], color=style.get('edge_fixed_color', 'blue'), figsize=(6, 6)))
        conditions.setdefault(condition, []).append((atlas, value))

    for condition, entries in conditions.items():
        labels = [atlas for atlas, _ in entries]
        folder = comparison_folder(entries[0][1])
        # Plotting must not perform undocumented inferential tests.
        common = dict(title=condition, xticklabel=labels, show_sub_label=True,
                      annotate_stats=False, figsize=(max(7, len(labels)*2), 4))
        add(folder + '/relative_effect.png', 'box',
            [np.asarray(v['relative_effect'], dtype=float) for _, v in entries],
            dict(ctx.viz_style('stat_box', 'stat_rel_effect_box'), **common, ylabel='1 - other / self'))
        group_style = dict(ctx.viz_style('stat_box', 'stat_group_indiv_box'), **common,
                           group_labels=['Other mice', 'Self'], ylabel='Pattern similarity (Fisher z)')
        group_style.update(show_sub_label=False, show_legend=True)
        add(folder + '/group_individual.png', 'box',
            [[v['other_z'] for _, v in entries], [v['self_z'] for _, v in entries]], group_style)

    for atlas in ['_'.join(a) for a in ctx.atlases]:
        for mode in ('observation_start', 'adjacent', 'windows'):
            entries = [v for key, v in complete.items()
                       if key.startswith(atlas + '/' + mode + ':') or key == atlas + '/windows' and mode == 'windows']
            if not entries:
                continue
            entries.sort(key=lambda v: (ctx.weeks.index(v['second'][0]), ctx.weeks.index(v['first'][0])))
            style = dict(ctx.viz_style('stat_mean_error'), title=f'{atlas}: {mode}',
                         xlabel='Observation comparison', ylabel='Pattern similarity (Fisher z)',
                         labels=['Other mice', 'Self'],
                         xticklabels=[','.join(v['first']) + ' → ' + ','.join(v['second']) for v in entries],
                         xtick_rotation=45, figsize=(max(9, len(entries)*.65), 5))
            style.update(legend_loc='upper center', legend_box=(.5, -.42))
            if len(entries) < 2 and style.get('diff_mode') == 'gradient':
                style.update(diff_mode='relative', diff_ylabel='Relative effect')
            add(f'report/figures/statistics/{atlas}/{mode}_group_individual.png', 'mean_error',
                ([v['other_z'] for v in entries], [v['self_z'] for v in entries]), style)
        pairs = [v for key, v in complete.items() if key.startswith(atlas + '/all_pairs:')]
        if pairs:
            for metric in ('accuracy', 'advantage_z'):
                matrix = np.full((len(ctx.weeks), len(ctx.weeks)), np.nan)
                for value in pairs:
                    i, j = ctx.weeks.index(value['first'][0]), ctx.weeks.index(value['second'][0])
                    matrix[j, i] = float(np.mean(value[metric]))
                add(f'report/figures/statistics/{atlas}/all_pairs_{metric}.png', 'edge_matrix',
                    (matrix, ctx.weeks), dict(title=f'{atlas}: {metric}', cmap='viridis', figsize=(8, 7),
                                             xlabel='First observation week', ylabel='Second observation week'))
        paired = {key: value for key, value in reliability.items()
                  if key.startswith(atlas + '/') and value.get('status') == 'complete'}
        if paired:
            add(f'report/figures/reliability/{atlas}/within_week.png', 'reliability', paired,
                dict(title=f'{atlas}: within-week identification', figsize=(10, 4)))
    structure = structure_metrics(ctx, results)
    if len(structure):
        groups = list(structure.group.unique())
        for metric in ('DD', 'ND'):
            arrays = [structure.loc[structure.group == group, metric].dropna().to_numpy() for group in groups]
            add(f'report/figures/statistics/structure/{metric}_box.png', 'box', arrays,
                dict(ctx.viz_style('stat_box', f'supplement_{metric.lower()}_box'),
                     title=f'{metric}: unique observation pairs (descriptive)', xticklabel=groups,
                     show_sub_label=True, annotate_stats=False, figsize=(9, 5)))
        data = [dict(x=structure.loc[structure.group == group, 'DD'].to_numpy(),
                     y=structure.loc[structure.group == group, 'ND'].to_numpy(), label=group) for group in groups]
        add('report/figures/statistics/structure/structure_space.png', 'scatter', data,
            dict(ctx.viz_style('supplement_structure_space'), xlabel='DD', ylabel='ND',
                 title='Structure space: unique observation pairs (descriptive)', figsize=(8, 6)))
    return jobs


def _draw(ctx, job):
    style = dict(job.style)
    figsize = style.pop('figsize', (6, 5))
    dpi = style.pop('dpi', ctx.figure_dpi())
    if job.kind == 'trajectory':
        from srcs.visualization.reports import plot_trajectories
        return plot_trajectories(job.data)
    if job.kind == 'circular':
        from srcs.visualization.circular_graph import circularGraph
        matrix, labels = job.data
        title = style.pop('title')
        if len(labels) % 2 or not all(str(x).endswith(('_L', '_R')) for x in labels):
            # Synthetic and non-hemispheric atlas labels need a general layout.
            return _generic_circle(matrix, labels, title, style, figsize)
        graph = circularGraph(matrix, [str(l).split('-', 1)[-1] for l in labels],
                              node_color_spec=labels, **style)
        graph.fig.set_size_inches(*figsize)
        graph.ax.set_aspect('equal')
        graph.ax.set_title(title, pad=25, fontsize=10)
        return graph.fig
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    plotter = PlotSeries()
    if job.kind in ('similarity', 'edge_matrix', 'edge_proportion'):
        matrix, labels = job.data
        if job.kind == 'edge_proportion':
            from matplotlib.colors import LinearSegmentedColormap
            matrix, labels = edge_proportions(matrix, labels, style['threshold'])
            style['cmap'] = LinearSegmentedColormap.from_list('edge_proportion', ['white', style['color']])
        matrix = np.array(matrix, dtype=float, copy=True)
        triangular = job.kind in ('edge_matrix', 'edge_proportion')
        if triangular:
            matrix[np.triu_indices(len(matrix), 0 if job.kind == 'edge_matrix' else 1)] = np.nan
        cmap = plt.get_cmap(style.get('cmap', 'inferno')).copy()
        cmap.set_bad('white')
        # Current values determine the color scale; legacy fixed ranges clip new results.
        artist = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, interpolation='none',
                           vmin=style.get('vmin'), vmax=style.get('vmax'))
        if style.get('cbar', True):
            bar = fig.colorbar(artist, ax=ax, shrink=.8)
            if job.kind == 'edge_proportion':
                bar.set_label('Edge proportion')
        fontsize = min(10, 300 / max(len(labels), 1))
        ax.set_xticks(np.arange(len(labels)), labels, rotation=90, fontsize=fontsize)
        ax.set_yticks(np.arange(len(labels)), labels, fontsize=fontsize)
        ax.set(title=style.get('title', ''), xlabel=style.get('xlabel', ''), ylabel=style.get('ylabel', ''))
        if triangular:
            lower_triangle_grid(ax, len(labels))
    elif job.kind == 'reliability':
        keys = list(job.data)
        x = np.arange(len(keys))
        ax.plot(x, [job.data[k]['accuracy'] for k in keys], 'o-', label='Bidirectional accuracy')
        ax.set_xticks(x, [k.split('/', 1)[1] for k in keys], rotation=45, ha='right')
        ax.set(title=style['title'], ylabel='Accuracy', ylim=(-.05, 1.05))
        ax.legend()
    else:
        if job.kind == 'hist':
            # Observed and null accuracies often coincide; label them in axes space.
            show_text = style.get('show_text', True)
            plotter._plot_hist(ax, job.data, dict(style, show_text=False))
            if show_text:
                ax.text(.97, .94, style.get('arrow_text_format', 'real={:.4f}').format(style['arrow_x']),
                        transform=ax.transAxes, ha='right', va='top', color=style.get('arrow_color', 'red'))
                ax.text(.97, .85, style.get('mean_text_format', 'shuffle={:.4f}').format(np.mean(job.data)),
                        transform=ax.transAxes, ha='right', va='top', color=style.get('mean_arrow_color', 'gray'))
        else:
            getattr(plotter, '_plot_' + job.kind)(ax, job.data, style)
    fig.tight_layout()
    return fig


def _generic_circle(matrix, labels, title, style, figsize):
    from matplotlib.patches import FancyArrowPatch
    fig, ax = plt.subplots(figsize=figsize)
    angle = np.linspace(0, 2*np.pi, len(labels), endpoint=False)
    xy = np.column_stack((np.cos(angle), np.sin(angle)))
    for i, j in zip(*np.triu_indices(len(labels), 1)):
        if matrix[i, j] != 0 and abs(matrix[i, j]) >= style.get('threshold', 0):
            ax.add_patch(FancyArrowPatch(xy[i], xy[j], connectionstyle='arc3,rad=.2',
                                        arrowstyle='-', color=style.get('edge_fixed_color', 'blue')))
    ax.scatter(*xy.T, zorder=3)
    for point, label in zip(xy, labels):
        ax.text(*(point*1.15), label, ha='center', va='center', fontsize=8)
    ax.set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), aspect='equal', title=title)
    ax.axis('off')
    return fig


def _render_file(config, job):
    from srcs.utils.config import RunContext
    ctx = RunContext(None, config, logging.getLogger('project_logger'))
    options = ctx.render_options()
    path = Path(ctx.output_path) / job.path
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        if job.kind == 'fc':
            from srcs.analysis.connectivity import CorrelationCache
            CorrelationCache.plot(SimpleNamespace(ctx=ctx), job.data, path.stem, resume=True)
        else:
            fig = _draw(ctx, job)
            if options['save']:
                _atomic(path, lambda f: fig.savefig(f, format='png', dpi=ctx.render_dpi(),
                        bbox_inches='tight', pil_kwargs={'compress_level': 1}), binary=True)
            if options['show']:
                plt.show()
        if options['save']:
            with Image.open(path) as img:
                img.verify()
            return hashlib.sha256(path.read_bytes()).hexdigest()
    finally:
        plt.close('all')


def _worker_init():
    plt.switch_backend('Agg')


def render_jobs(ctx, jobs, index_path, workers=1):
    """Per-file verified resume, including partial recovery after a failed report."""
    root = Path(ctx.output_path)
    options = ctx.render_options()
    previous = {}
    if Path(index_path).exists():
        try:
            previous = {row['path']: row for row in json.loads(Path(index_path).read_text(encoding='utf-8'))['figures']}
        except (ValueError, KeyError):
            pass
    code_hash = hashlib.sha256(b''.join(p.read_bytes() for p in sorted(Path(__file__).parent.glob('*.py')))).hexdigest()
    index = dict(status='running', config_tag=ctx.analysis_tag, figures=[],
                 interpretation='Descriptive figures. No implicit hypothesis tests or multiple-comparison correction.',
                 save=options['save'], show=options['show'])
    save_json(index_path, index)
    if not options['save'] and not options['show']:
        index['status'] = 'disabled'
        save_json(index_path, index)
        return
    pool = None
    try:
        pending = []
        for job in progress_items(jobs, phase='checking figure cache'):
            path = root / job.path
            identity = hashlib.sha256(pickle.dumps((job.kind, job.data, job.style, ctx.render_dpi(), code_hash))).hexdigest()
            old = previous.get(job.path, {})
            reusable = (options['save'] and not options['show'] and old.get('identity') == identity and path.is_file()
                        and old.get('sha256') == hashlib.sha256(path.read_bytes()).hexdigest())
            if reusable:
                index['figures'].append(old)
            else:
                pending.append((job, identity))
        completed = len(jobs) - len(pending)
        ctx.logger.info('Report figures: %s/%s verified; %s to render', completed, len(jobs), len(pending))
        emit_progress(completed, len(jobs), 'figures')

        def record(job, identity, sha256):
            nonlocal completed
            if options['save']:
                index['figures'].append(dict(path=job.path, kind=job.kind, identity=identity, sha256=sha256))
            completed += 1
            emit_progress(completed, len(jobs), 'figures')
            if completed % 20 == 0:
                save_json(index_path, index)
                ctx.logger.info('Report figures: %s/%s', completed, len(jobs))

        if workers > 1 and not options['show'] and pending:
            pool = ProcessPoolExecutor(max_workers=workers, initializer=_worker_init)
            futures = {pool.submit(_render_file, ctx.config, job): (job, identity) for job, identity in pending}
            for future in as_completed(futures):
                job, identity = futures[future]
                record(job, identity, future.result())
        else:
            for job, identity in pending:
                record(job, identity, _render_file(ctx.config, job))
    except BaseException:
        index['status'] = 'incomplete'
        save_json(index_path, index)
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    order = {job.path: i for i, job in enumerate(jobs)}
    index['figures'].sort(key=lambda row: order[row['path']])
    index['status'] = 'complete'
    save_json(index_path, index)
