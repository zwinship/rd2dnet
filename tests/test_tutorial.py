"""Execute the tutorial notebook end to end (run with ``pytest -m notebook``)."""

from pathlib import Path

import pytest

NOTEBOOK = Path(__file__).resolve().parents[1] / "examples" / "tutorial.ipynb"


@pytest.mark.notebook
def test_tutorial_runs():
    nbformat = pytest.importorskip("nbformat")
    nbclient = pytest.importorskip("nbclient")
    nb = nbformat.read(NOTEBOOK, as_version=4)
    nbclient.NotebookClient(nb, timeout=600, kernel_name="python3",
                            resources={"metadata": {"path": str(NOTEBOOK.parent)}}).execute()
