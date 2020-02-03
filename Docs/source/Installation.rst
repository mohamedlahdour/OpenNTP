Quick Install Guide
=============

This quick install guide outlines the basic steps needed to install OpenNTP on your computer.

Installing on Linux
*******************

1. If you are using Ubuntu 18.10 or newer, open a terminal in a GNU/Linux box then install gfortran with the following commands:
    .. code-block:: python

        sudo apt-get update
        sudo apt-get install gfortran

2. You need to install numpy (F2PY) and matplotlib library to run the package OpenNTP:
    .. code-block:: python

        sudo apt install python-pip
        sudo apt-get install python-numpy
        pip3 install numpy matplotlib

3. You need to install PyQt5 on Ubuntu with python to run the GUI:
    .. code-block:: python

        pip install --user pyqt5  
        sudo apt-get install python-pyqt5  
        sudo apt-get install pyqt5-dev-tools
        sudo apt-get install qttools5-dev-tools

4. Install the OpenNTP package
    .. code-block:: python

        git clone  https://github.com/mohamedlahdour/OpenNTP.git

5. Import the *OpenNTP* and run the package in the following way:
    .. code-block:: python

         cd OpenNTP
         $ python main.py
