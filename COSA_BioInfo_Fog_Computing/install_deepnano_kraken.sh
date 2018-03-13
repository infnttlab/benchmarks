#!/bin/bash

export LC_ALL="en_US.UTF-8"
export LC_CTYPE="en_US.UTF-8"

mkdir fog_computing
cd fog_computing/
    git clone https://elenacorni@bitbucket.org/vboza/deepnano.git
    
    sudo apt-get install python-pip
    sudo -H python -m pip install Theano
    sudo -H python -m pip install h5py
    sudo -H python -m pip install python-dateutil
    
    cd deepnano/
        g++ -O2 -std=gnu++0x align_2d.cc -o align_2d
    cd ..
    
    wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz
    tar -zxvf jellyfish-1.1.11.tar.gz

    cd jellyfish-1.1.11
        ./configure
        make
        make install
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
    cd ..
    
    wget https://ccb.jhu.edu/software/kraken/dl/kraken-0.10.5-beta.tgz
    tar -xvzf kraken-0.10.5-beta.tgz

    cd kraken-0.10.5-beta
        ./install_kraken.sh /root/fog_computing/kraken-0.10.5-beta
    cd ..
    
    wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
    tar -xvzf minikraken.tgz
    
    rm minikraken.tgz kraken-0.10.5-beta.tgz jellyfish-1.1.11.tar.gz
cd ..
echo "Done."
