"""Functions to WEAS visualisation."""

from __future__ import annotations

from pathlib import Path


def generate_weas_html(filename: str | Path) -> str:
    """
    Generate HTML for WEAS.

    Parameters
    ----------
    filename
        Path of structure file.

    Returns
    -------
    str
        HTML for WEAS to visualise structure.
    """
    return f"""
    <!doctype html>
    <html lang="en">
    <body>
        <div id="viewer" style="position: relative; width: 100%; height: 500px"></div>

        <script type="module">

        async function fetchFile(filename) {{
            const response = await fetch(`${{filename}}`);
            if (!response.ok) {{
            throw new Error(`Failed to load file for structure: ${{filename}}`);
            }}
            return await response.text();
        }}

        import {{ WEAS, parseXYZ, parseCIF, parseCube, parseXSF }} from 'https://unpkg.com/weas/dist/index.mjs';
        const domElement = document.getElementById("viewer");

        // hide the buttons
        const guiConfig = {{
            buttons: {{
                enabled: false,
            }},
        }};
        const editor = new WEAS({{ domElement, viewerConfig: {{ _modelStyle: 1 }}, guiConfig}});

        let structureData;
        const filename = "{str(filename)}";
        console.log("filename: ", filename);
        structureData = await fetchFile(filename);
        console.log("structureData: ", structureData);

        if (filename.endsWith(".xyz") || filename.endsWith(".extxyz")) {{

            const atoms = parseXYZ(structureData);
            editor.avr.atoms = atoms;
            editor.avr.modelStyle = 1;

        }} else if (filename.endsWith(".cif")) {{

            const atoms = parseCIF(structureData);
            editor.avr.atoms = atoms;
            editor.avr.showBondedAtoms = true;
            editor.avr.colorType = "VESTA";
            editor.avr.boundary = [[-0.01, 1.01], [-0.01, 1.01], [-0.01, 1.01]];
            editor.avr.modelStyle = 2;

        }} else if (filename.endsWith(".cube")) {{

            const data = parseCube(structureData);
            editor.avr.atoms = data.atoms;
            editor.avr.volumetricData = data.volumetricData;
            editor.avr.isosurfaceManager.fromSettings({{
            positive: {{ isovalue: 0.0002, mode: 1, step_size: 1 }},
            negative: {{ isovalue: -0.0002, color: "#ff0000", mode: 1 }},
            }});
            editor.avr.isosurfaceManager.drawIsosurfaces();

        }} else if (filename.endsWith(".xsf")) {{

            const data = parseXSF(structureData);
            editor.avr.atoms = data.atoms;
            editor.avr.volumetricData = data.volumetricData;
            editor.avr.isosurfaceManager.fromSettings({{
            positive: {{ isovalue: 0.15, mode: 1 }},
            negative: {{ isovalue: -0.15, color: "#ff0000", mode: 1 }},
            }});
            editor.avr.isosurfaceManager.drawIsosurfaces();

        }} else {{
            document.getElementById("viewer").innerText = "Unsupported file format.";
        }}

        editor.render();

        </script>
    </body>
    </html>
    """  # noqa: E501
