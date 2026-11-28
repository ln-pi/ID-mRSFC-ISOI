"""Publish unified figures as HTML, PDF and a checksummed manifest."""
import base64
import html
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from srcs.core.artifacts import artifact_ref, write_json


def publish_report(ctx, figures, inputs, skipped):
    folder = Path(ctx.output_dir("report"))
    entries, body = [], []
    with PdfPages(folder/"report.pdf") as pdf:
        for path in figures:
            path = Path(path)
            entries.append(artifact_ref(path, folder))
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            title = html.escape(path.stem)
            body.append(f'<figure><img alt="{title}" src="data:image/png;base64,{encoded}"><figcaption>{title}</figcaption></figure>')
            fig,ax=plt.subplots(figsize=(11.7,8.3))
            ax.imshow(plt.imread(path))
            ax.axis("off")
            pdf.savefig(fig,bbox_inches="tight")
            plt.close(fig)
        if not figures:
            fig,ax=plt.subplots()
            ax.text(.1,.5,"No figures available")
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)
    (folder/"report.html").write_text('<!doctype html><meta charset="utf-8"><title>FC analysis</title>'
        '<style>body{font:16px system-ui;max-width:1100px;margin:auto}img{max-width:100%}</style>'
        '<h1>FC 分析报告</h1>'+''.join(body),encoding="utf-8")
    write_json(folder/"manifest.json",dict(figures=entries,inputs=inputs,skipped=skipped))
