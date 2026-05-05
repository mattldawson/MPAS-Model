#!/bin/bash
set -euo pipefail

# Install MUSICA
cd /musica-src
cmake --install build-container 2>&1 | tail -1
ldconfig

cd /mpas

# No-op external checkout scripts
printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/tools/manage_externals/checkout_externals
chmod +x src/core_atmosphere/tools/manage_externals/checkout_externals
printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/physics/checkout_data_files.sh
chmod +x src/core_atmosphere/physics/checkout_data_files.sh

# Clean
echo y | make clean CORE=atmosphere 2>/dev/null || true

# Build
JLEVEL=$(( $(nproc) > 8 ? 8 : $(nproc) ))
make -j${JLEVEL} gnu CORE=atmosphere USE_PIO2=false \
    MPAS_EXTERNAL_LIBS="$(pkg-config --libs musica-fortran) -lstdc++" \
    MPAS_EXTERNAL_INCLUDES="$(pkg-config --cflags musica-fortran)" 2>&1 | tail -5

echo "=== BUILD RESULT ==="
ls -la atmosphere_model
echo "DONE"
