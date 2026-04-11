sudo apt update
# install build-essential if not already installed
if ! dpkg -s build-essential &> /dev/null
then
    echo "build-essential could not be found, installing..."
    sudo apt install -y build-essential
else
    echo "build-essential is already installed"
fi
# install fpc only if it not already installed
if ! command -v fpc &> /dev/null
then
    echo "fpc could not be found, installing..."
    sudo apt install -y fpc
else
    echo "fpc is already installed"
fi
# install git, wget, tar if not already installed
if ! command -v git &> /dev/null
then
    echo "git could not be found, installing..."
    sudo apt install -y git
else
    echo "git is already installed"
fi
if ! command -v wget &> /dev/null
then
    echo "wget could not be found, installing..."
    sudo apt install -y wget
else
    echo "wget is already installed"
fi
if ! command -v tar &> /dev/null
then
    echo "tar could not be found, installing..."
    sudo apt install -y tar
else
    echo "tar is already installed" 
fi
cd ..
pwd
# clone core-math repository if not already cloned
if [ ! -d "core-math" ]; then
    git clone https://gitlab.inria.fr/core-math/core-math.git
fi
pwd
