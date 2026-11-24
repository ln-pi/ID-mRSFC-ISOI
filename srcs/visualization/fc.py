"""Optional per-acquisition FC diagnostic rendering."""
import hashlib
import json
from pathlib import Path
from PIL import Image
from srcs.visualization.plot_series import save_fc_heatmap
from srcs.utils.storage import save_json
from srcs.visualization.scales import resolve


def render_fc(ctx, frame, stem, resume=True):
    options, style = ctx.render_options(), ctx.viz_style("fc_heatmap")
    viz = ctx.viz()
    scale = resolve(frame.to_numpy(), viz, "fc_heatmap", bounds=(-1, 1))
    if scale is None:
        return
    for key in ("center", "vmin", "vmax"):
        style.pop(key, None)
    style.update(cbar_extend=scale.extend)
    if not options.get("save", True) and not options.get("show", False):
        return
    folder = Path(ctx.output_dir("Correlation", "Figures"))
    png, marker = folder / (stem + ".png"), folder / (stem + ".plot.json")
    identity = {"values": hashlib.sha256(frame.to_numpy().tobytes()).hexdigest(),
                "labels": list(frame.index), "options": options, "style": style,
                "scale": scale.record(), "special_colors": {
                    key: viz[key] for key in ("missing_color", "under_color", "over_color")}}
    complete = False
    if resume and png.exists() and marker.exists():
        try:
            complete = json.loads(marker.read_text(encoding="utf-8")) == identity
            with Image.open(png) as img:
                img.verify()
        except Exception:
            complete = False
    if complete and not options.get("show", False):
        return
    save_fc_heatmap(frame.to_numpy(), names=list(frame.index), name=stem,
                    save_path=str(folder), config={**style, **scale.image}, **options)
    if options.get("save", True):
        with Image.open(png) as img:
            img.verify()
        save_json(marker, identity)
