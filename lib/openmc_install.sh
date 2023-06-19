#Begin.

# Install dependancies.
cd ~
sudo apt update
sudo apt install git cmake g++ mpich gfortran libhdf5-dev libpng-dev libeigen3-dev libmpich-dev libnetcdf-dev libtbb-dev libgles2-mesa-dev
pip install -r openmc_requirements.txt

#Install openmc.
git clone https://github.com/openmc-dev/openmc.git
#git clone https://gitclone.com/github.com/openmc-dev/openmc.git
cd openmc
pip install . && cmake .
sudo make install
cd ~

#End.
