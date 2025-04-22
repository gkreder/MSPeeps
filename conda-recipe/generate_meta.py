#!/usr/bin/env python
"""
Script to generate meta.yaml from pyproject.toml.
Run this script from the project root to update conda-recipe/meta.yaml.
"""
import tomli
import os
import sys

# Read pyproject.toml
with open("pyproject.toml", "rb") as f:
    pyproject = tomli.load(f)

# Extract metadata
project = pyproject.get("project", {})
name = project.get("name", "mspeeps")
version = project.get("version", "0.1.0")
description = project.get("description", "Mass spectrometry data extraction and conversion")
license = project.get("license", {}).get("text", "MIT")
python_requires = project.get("requires-python", ">=3.10")
dependencies = project.get("dependencies", [])

# Format dependencies for conda
conda_deps = []
for dep in dependencies:
    # Convert from PEP 508 format to conda format if needed
    if ";" in dep:  # Skip environment markers for now
        continue
    
    # Handle pip vs conda package names (potential differences)
    if "rdkit" in dep.lower():
        conda_deps.append("rdkit")
        continue
    
    # Handle version specifiers
    if ">=" in dep and "<" in dep:
        name_part = dep.split(">=")[0].strip()
        min_ver = dep.split(">=")[1].split("<")[0].strip()
        max_ver = dep.split("<")[1].strip()
        conda_deps.append(f"{name_part} >={min_ver},<{max_ver}")
    elif ">=" in dep:
        parts = dep.split(">=")
        name_part = parts[0].strip()
        version = parts[1].strip().split(",")[0].strip()
        conda_deps.append(f"{name_part} >={version}")
    elif "==" in dep:
        parts = dep.split("==")
        name_part = parts[0].strip()
        version = parts[1].strip().split(",")[0].strip()
        conda_deps.append(f"{name_part} =={version}")
    else:
        conda_deps.append(dep)

# Create meta.yaml content
meta_yaml = f"""package:
  name: {name}
  version: "{version}"

source:
  path: ..

build:
  noarch: python
  script: python -m pip install .
  entry_points:
    - mspeeps = mspeeps.cli:main

requirements:
  host:
    - python {python_requires}
    - pip
    - setuptools
  run:
    - python {python_requires}
"""

# Add dependencies
for dep in conda_deps:
    meta_yaml += f"    - {dep}\n"

meta_yaml += f"""
about:
  home: https://github.com/gkreder/mspeeps
  summary: "{description}"
  license: {license}
  description: |
    MSPeeps is a Python package for extracting mass spectrometry spectra
    from mzML files and converting them to MSP format. It provides both
    a command-line interface and a Python API.
  dev_url: https://github.com/gkreder/mspeeps
  doc_url: https://github.com/gkreder/mspeeps
  doc_source_url: https://github.com/gkreder/mspeeps/blob/main/README.md
"""

# Write to meta.yaml
with open("conda-recipe/meta.yaml", "w") as f:
    f.write(meta_yaml)

print("Generated conda-recipe/meta.yaml from pyproject.toml") 