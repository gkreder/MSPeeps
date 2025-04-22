#!/usr/bin/env python
"""
Script to generate meta.yaml from pyproject.toml.
Run this script from the project root to update conda-recipe/meta.yaml.
"""
import tomli
import os
import sys
import re

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
    # Skip environment markers
    if ";" in dep:
        continue
    
    # Extract name and version constraints using regex
    name_match = re.match(r"([a-zA-Z0-9_\-]+)", dep)
    name_part = name_match.group(1) if name_match else dep
    
    # Handle rdkit specifically
    if name_part.lower() == "rdkit":
        conda_deps.append("rdkit")
        continue
    
    # Special case for pandas or any package with complex version range
    if name_part.lower() == "pandas":
        # Based on pyproject.toml, specify a fixed format that conda accepts
        conda_deps.append(f"{name_part} >=2.2.3")
        continue
    
    # Handle version specifiers more robustly for other packages
    version_part = dep[len(name_part):].strip()
    
    if version_part:
        # Format version constraints properly for conda
        if ">=" in version_part and "<" in version_part:
            # Extract min and max versions
            min_ver_match = re.search(r">=\s*([0-9.]+)", version_part)
            max_ver_match = re.search(r"<\s*([0-9.]+)", version_part)
            
            min_ver = min_ver_match.group(1) if min_ver_match else None
            max_ver = max_ver_match.group(1) if max_ver_match else None
            
            if min_ver and max_ver:
                conda_deps.append(f"{name_part} >={min_ver},<{max_ver}")
            elif min_ver:
                conda_deps.append(f"{name_part} >={min_ver}")
            elif max_ver:
                conda_deps.append(f"{name_part} <{max_ver}")
            else:
                conda_deps.append(name_part)
        else:
            # Simple version constraint (like ==, >=, etc.)
            version_part = re.sub(r',\s*', ',', version_part)  # Remove spaces after commas
            conda_deps.append(f"{name_part} {version_part}")
    else:
        conda_deps.append(name_part)

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