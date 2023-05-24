export WORKDIR=$(pwd)

function load_modules() {
    # module purge
    # module load slurm
    # module load dot
    # module load turbovnc/2.0.1
    # module load vgl/2.5.1/64
    # module load singularity/current
    # module load rhel7/global
    # module load cmake/latest
    # module load gcc/7
    # module load openmpi-3.1.3-gcc-7.2.0-b5ihosk
    # module load python/3.7
    # module load hdf5/1.12.1
    # export LIB_PATH=""

    export CC=mpicc
    export CXX=mpic++
    export FF=mpif90

    source ~/.python_venv/ML/bin/activate
}

function install_libxc(){
    git clone https://gitlab.com/libxc/libxc.git
    cd libxc

    mkdir build
    mkdir bin
    cd build

    cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FF -DENABLE_FORTRAN=ON -DCMAKE_INSTALL_PREFIX=$WORKDIR/libxc/bin ..
    make -j5
    # make test
    make install

    cd $WORKDIR/libxc

    python3 setup.py install

    # export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORKDIR/libxc/bin/lib
    # export LD_LIBRARY_PATH=$WORKDIR/libxc/bin/lib

    # export PATH=$PATH:$WORKDIR/libxc/bin

    cd $WORKDIR
}


function install_QE(){

    git clone https://gitlab.com/QEF/q-e.git
    cd q-e
    mkdir build
    mkdir bin

    cd build
    
    # cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DQE_ENABLE_LIBXC=ON ..
    # cmake -DCMAKE_INSTALL_PREFIX=$WORKDIR/q-e/bin -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FF -DQE_LAPACK_INTERNAL=ON -DQE_FFTW_VENDOR=Internal -DQE_ENABLE_LIBXC=ON ..

    export Libxc_DIR=$WORKDIR/libxc/bin
    cmake -DCMAKE_INSTALL_PREFIX=$WORKDIR/q-e/bin -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FF -DQE_LAPACK_INTERNAL=ON -DQE_FFTW_VENDOR=Internal -DQE_ENABLE_LIBXC=ON ..
    
    # CC=$CC CXX=$CXX FC=$FC ./configure --with-libxc --with-libxc-prefix='/Users/akgray/Documents/DFT/libxc/bin' -–with-libxc-include='/Users/akgray/Documents/DFT/libxc/bin/include'
    CC=$CC CXX=$CXX FC=$FC ./configure --with-libxc --with-libxc-prefix="$WORKDIR/Users/akgray/Documents/DFT/libxc/bin" --prefix="/Users/akgray/Documents/DFT/q-e"
    make -j5 all
    make install

    cd $WORKDIR
}


function install_ASE(){     
    git clone -b 3.22.1 https://gitlab.com/ase/ase.git
    cd ase

    python3 setup.py install
    # ase test

    cd $WORKDIR

    git clone https://github.com/lmmentel/ase-espresso.git
    cd ase-espresso
    
    pip3 install .

    export PATH="$WORKDIR/q-e/bin:"$PATH       # add these to bashrc
    export ESP_PSP_PATH="$WORKDIR/q-e/pseudo"  # add these to bashrc

    cd $WORKDIR
}