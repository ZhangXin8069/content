home='/content/drive/MyDrive/'
cd $home
pip install --upgrade pip
sudo apt update
sudp apt install git cmake g++ mpich gfortran libhdf5-dev libpng-dev libeigen3-dev libmpich-dev libnetcdf-dev libtbb-dev libgles2-mesa-dev
pip install -r openmc_requirements.txt
cd openmc
sudo make install
pip install . && cmake .
cd $home
python colab_prepare_openmc.py
