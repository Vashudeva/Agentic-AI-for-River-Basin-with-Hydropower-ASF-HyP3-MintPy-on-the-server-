#!/usr/bin/env bash
# Create the kgproc environment in Linux / WSL2.
#   bash setup_env.sh [ENV_PREFIX]        (default: ~/kali-insar-env)
# Needs internet access to conda-forge and GitHub. Takes 10-20 minutes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="${1:-$HOME/kali-insar-env}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

if ! command -v micromamba >/dev/null 2>&1; then
  echo "Installing micromamba into $HOME/bin"
  mkdir -p "$HOME/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME" bin/micromamba
  export PATH="$HOME/bin:$PATH"
fi

echo "Creating environment $ENV_PREFIX"
micromamba create -y -p "$ENV_PREFIX" -f "$HERE/environment.yml"
PY="$ENV_PREFIX/bin/python"

echo "Installing MiaplPy (GitHub, not on conda-forge)"
"$PY" -m pip install --no-deps "git+https://github.com/insarlab/MiaplPy.git"

echo "Linking the snaphu executable (MiaplPy calls 'snaphu')"
SNAPHU="$("$PY" -c 'import snaphu, os; print(os.path.join(os.path.dirname(snaphu.__file__), "snaphu"))')"
ln -sf "$SNAPHU" "$ENV_PREFIX/bin/snaphu"

echo "Activation hook: ISCE2 topsStack, ISCE components and MiaplPy scripts"
mkdir -p "$ENV_PREFIX/etc/conda/activate.d"
SITE="$("$PY" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
cat > "$ENV_PREFIX/etc/conda/activate.d/kgproc.sh" <<EOF
export PATH="\$CONDA_PREFIX/share/isce2/topsStack:$SITE/isce/applications:\$PATH"
export PYTHONPATH="\$CONDA_PREFIX/share/isce2:$SITE/isce/components:$SITE/miaplpy\${PYTHONPATH:+:\$PYTHONPATH}"
EOF

echo
(cd "$HERE" && "$PY" -m kgproc check -c config.toml) || true
echo
echo "Done. Run the pipeline with:"
echo "  cd $HERE && $PY -m kgproc all"
