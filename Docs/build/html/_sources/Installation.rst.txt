Quick Install Guide
=============

This quick install guide outlines the basic steps needed to install **OpenNTP** on your computer.

Installing on Linux
*******************

1. If you are using Ubuntu 20.04, open a terminal in a GNU/Linux box then install gfortran with the following commands:
    .. code-block:: python

        sudo apt-get install gfortran

2. You need to install numpy (F2PY3) and matplotlib library to run the package OpenNTP:
    .. code-block:: python

        sudo apt-get install python3-numpy
        sudo apt-get install python3-matplotlib

3. You need to install PyQt5 on Ubuntu with python to run the GUI:
    .. code-block:: python

        sudo apt-get install python3-pyqt5 

4. Install the **OpenNTP** package
    .. code-block:: python

        git clone  https://github.com/mohamedlahdour/OpenNTP.git

5. Import the **OpenNTP** and run the package in the following way:
    .. code-block:: python

        cd OpenNTP
        $ python3 gui.py
