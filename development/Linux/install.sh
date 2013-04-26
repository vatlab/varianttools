#!/usr/bin/env bash
set -e
is_relative() {
    local path="$1"
    shift

    [ "${path:0:1}" != "/" ]
    return
}
install () {
    mkdir -p $2/bin
    rm -rf $2/lib/variant_tools &> /dev/null
    mkdir -p $2/lib/variant_tools
    cp -a $1/* $2/lib/variant_tools
    for cmd in vtools vtools_report; do
	rm -rf $2/bin/$cmd &> /dev/null
	ln -s $2/lib/variant_tools/$cmd $2/bin/$cmd
    done
    echo -e "Libraries are installed to $2/lib\nBinary files are installed to $2/bin\n"
}
main () {
    local fullpath=""
    echo "Enter installation directory for variant tools & variant association tools: "
    printf "\t [/usr/local]  "
    read fullpath
    if [ -z $fullpath ]; then
	install $1 "/usr/local"
    else
	eval fullpath=$fullpath
	if is_relative $fullpath; then
	    fullpath=$PWD/$fullpath
	fi
	install $1 $fullpath
    fi
}
main $@
