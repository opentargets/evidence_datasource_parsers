#!/bin/bash

if ! command -v uv &>/dev/null; then
    echo "uv was not found, installing uv..."
    wget -qO- https://astral.sh/uv/install.sh | sh
    source "$HOME/.${SHELL##*/}rc"
fi

# Install the required Python version defined in pyproject.toml
PYTHON_VERSION=$(grep 'requires-python' pyproject.toml | awk -F'"' '{print $2}' | cut -d'>=' -f2 | cut -d'<' -f1)
echo "Installing Python version $PYTHON_VERSION..."
uv python install "$PYTHON_VERSION"

# Create a virtual environment if not already present
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv venv .venv
fi

# Activate the virtual environment and install dependencies
echo "Activating virtual environment..."
source .venv/bin/activate
echo "Installing dependencies with UV..."
uv sync --all-groups --frozen
echo "Installing project in editable mode..."
uv pip install -e .