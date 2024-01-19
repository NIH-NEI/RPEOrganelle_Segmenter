cd swigs\rpesegm
python setup.py clean build install_lib --force --install-dir ..\..
copy /Y rpesegm.py ..\..
copy /Y rpeutil.py ..\..
cd ..\..
