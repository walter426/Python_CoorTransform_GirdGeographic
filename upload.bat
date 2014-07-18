setlocal
set HOME=C:\Python27
#python setup.py register
python setup.py sdist bdist_wininst upload
endlocal
pause