from pathlib import Path

import subprocess

import julia


here = Path(__file__).parent.resolve()

def install() -> None:
    """Install package."""
    julia.install()

    subprocess.run(f"python3 -m julia.sysimage {here}/sys.so", shell=True)