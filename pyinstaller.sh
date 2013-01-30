for exe in vtools vtools_report; do
    python $@/pyinstaller.py -F $exe 
done
