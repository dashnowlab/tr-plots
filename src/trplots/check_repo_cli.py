import ast
import os
from pathlib import Path

# files to scan root (adjust if needed)
ROOT = Path(__file__).resolve().parents[2]  # project root -> .../tr-plots
SRC_DIR = ROOT / "src" / "trplots"

CONFIG_NAMES = {
    "VCF_PATH", "JSON_PATH", "OUTPUT_BASE", "OTHER_DATA", "DATA_PATH",
    "OUTPUT_DIR", "OUTPUT_HTML_DIR", "OUTPUT_PNG_DIR", "CSV_PATH",
    "ALLELE_LENGTH_PLOTS_OUTPUT"
}

def analyze_file(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore")
    try:
        tree = ast.parse(text)
    except Exception:
        return {"path": str(path), "parse_error": True}

    info = {
        "path": str(path),
        "uses_trplots_config_import": False,
        "references_config_names": sorted([n for n in CONFIG_NAMES if n in text]),
        "has_main_def": any(isinstance(n, ast.FunctionDef) and n.name == "main" for n in tree.body),
        "has_parse_args_def": any(isinstance(n, ast.FunctionDef) and n.name == "parse_args" for n in tree.body),
        "has_main_guard": "__name__" in text and '== "__main__"' in text,
        "imports_argparse": any(isinstance(node, (ast.Import, ast.ImportFrom)) and 
                                ((isinstance(node, ast.Import) and any(alias.name == "argparse" for alias in node.names)) or
                                 (isinstance(node, ast.ImportFrom) and node.module == "argparse"))
                                for node in ast.walk(tree)),
        "imports_plotly_io": any(isinstance(node, (ast.Import, ast.ImportFrom)) and 
                                ((isinstance(node, ast.Import) and any(alias.name in ("plotly.io","plotly.io as pio","pio") for alias in node.names)) or
                                 (isinstance(node, ast.ImportFrom) and node.module == "plotly.io"))
                                for node in ast.walk(tree)),
        "sets_pio_renderer": "pio.renderers.default" in text,
        "fig_show_with_renderer": "fig.show(" in text and "renderer=" in text
    }

    # check for `from trplots.config import` or `import trplots.config`
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("trplots.config"):
            info["uses_trplots_config_import"] = True
            break
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("trplots.config") or alias.name == "trplots.config":
                    info["uses_trplots_config_import"] = True
                    break

    return info

def main():
    if not SRC_DIR.exists():
        print(f"ERROR: expected source dir not found: {SRC_DIR}")
        return

    py_files = sorted(SRC_DIR.rglob("*.py"))
    results = [analyze_file(p) for p in py_files]

    # Print summary table
    hdr = ["file", "uses_config_import", "config_refs", "has_main", "has_parse_args",
           "has_main_guard", "imports_argparse", "sets_pio_renderer", "fig_show_renderer"]
    print("\t".join(hdr))
    for r in results:
        if r.get("parse_error"):
            print(f"{r['path']}\tPARSE_ERROR")
            continue
        row = [
            os.path.relpath(r["path"], ROOT),
            str(r["uses_trplots_config_import"]),
            ",".join(r["references_config_names"]) or "-",
            str(r["has_main_def"]),
            str(r["has_parse_args_def"]),
            str(r["has_main_guard"]),
            str(r["imports_argparse"]),
            str(r["sets_pio_renderer"]),
            str(r["fig_show_with_renderer"])
        ]
        print("\t".join(row))

if __name__ == "__main__":
    main()