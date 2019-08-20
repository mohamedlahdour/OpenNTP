Quick Install Guide
=============

This quick install guide outlines the basic steps needed to install OpenRSN on your computer.

Installing on Linux
*******************

1. If you are using Ubuntu 18.10 or newer, then you can easily install Python 3 with the following commands:
    .. code-block:: python

        sudo apt-get update
        sudo apt-get install python3

2. Open a terminal in a GNU/Linux box and install the following tools:
    .. code-block:: python


        sudo apt-get install gfortran

3. You need to install numpy and matplotlib library to run the package OpenRSN:
    .. code-block:: python

        sudo apt-get install python3-pip  
        sudo pip3 install numpy 
        pip3 install numpy matplotlib

4. You need to install PyQt5 on Ubuntu with python3 to run the GUI:
    .. code-block:: python

        pip3 install --user pyqt5  
        sudo apt-get install python3-pyqt5  
        sudo apt-get install pyqt5-dev-tools
        sudo apt-get install qttools5-dev-tools

5. Install the OpenRSN package
    .. code-block:: python

        git clone  https://github.com/mohamedlahdour/OpenRSN.git

6. Import the *OpenRSN* and run the package in the following way:
    .. code-block:: python

         cd OpenRSN
         $ python3 main.py
