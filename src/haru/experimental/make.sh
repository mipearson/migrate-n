#/bin/sh
gcc -o DocMaker DocMaker.cc PdfWMFLoader.cc -I.. -L.. -lharu -DDEBUG \
    -lz -lstdc++ -lpng -DDEBUG

if [ $? -ne 0 ]; then
    echo "ERROR!!"
    exit -1
fi
    
./DocMaker ./libharu.txt libharu.pdf > /dev/null
    
exit $?

