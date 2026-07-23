from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml


def load_yaml(path: Path) -> dict[str, Any]:
    """Carga una configuración YAML y valida su estructura básica."""
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict):
        raise ValueError(
            f"Configuración YAML inválida: {path}"
        )

    return data


def write_yaml(
    path: Path,
    data: dict[str, Any],
) -> None:
    """Escribe un YAML de forma determinista."""
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            data,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )


def require_sections(
    config: dict[str, Any],
    sections: list[str],
) -> None:
    """Comprueba que existan las secciones requeridas."""
    missing = [
        section
        for section in sections
        if section not in config
    ]

    if missing:
        raise KeyError(
            "Faltan secciones obligatorias en la configuración: "
            + ", ".join(missing)
        )


def select_sections(
    config: dict[str, Any],
    sections: list[str],
) -> dict[str, Any]:
    """Selecciona únicamente las secciones usadas por un módulo."""
    require_sections(
        config=config,
        sections=sections,
    )

    return {
        section: config[section]
        for section in sections
    }


def build_load_config(
    config: dict[str, Any],
) -> dict[str, Any]:
    return select_sections(
        config,
        [
            "project",
            "input",
            "output",
        ],
    )


def build_qc_filter_config(
    config: dict[str, Any],
) -> dict[str, Any]:
    return select_sections(
        config,
        [
            "project",
            "output",
            "qc",
        ],
    )


def build_qc_visualization_config(
    config: dict[str, Any],
) -> dict[str, Any]:
    return select_sections(
        config,
        [
            "project",
            "output",
        ],
    )


def build_doublet_config(
    config: dict[str, Any],
) -> dict[str, Any]:
    return select_sections(
        config,
        [
            "project",
            "output",
            "doublet",
        ],
    )


def build_qc_scores_config(
    config: dict[str, Any],
) -> dict[str, Any]:
    return select_sections(
        config,
        [
            "project",
            "output",
            "scoring",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Genera configuraciones runtime independientes "
            "para los módulos del pipeline scRNA-seq."
        )
    )

    parser.add_argument(
        "--config-file",
        required=True,
        help="Configuración completa del proyecto.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directorio donde se escribirán los YAML runtime.",
    )

    args = parser.parse_args()

    config_file = Path(
        args.config_file
    ).expanduser().resolve()

    output_dir = Path(
        args.output_dir
    ).expanduser().resolve()

    if not config_file.is_file():
        raise FileNotFoundError(
            f"No existe la configuración: {config_file}"
        )

    config = load_yaml(config_file)

    runtime_configs = {
        "load_config.yaml": build_load_config(config),
        "qc_filter_config.yaml": build_qc_filter_config(config),
        "qc_visualization_config.yaml":
            build_qc_visualization_config(config),
        "doublet_config.yaml": build_doublet_config(config),
        "qc_scores_config.yaml": build_qc_scores_config(config),
    }

    for filename, runtime_config in runtime_configs.items():
        output_path = output_dir / filename

        write_yaml(
            path=output_path,
            data=runtime_config,
        )

        print(output_path)


if __name__ == "__main__":
    main()
