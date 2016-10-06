display_usage(){
        echo -e "\nUsage: $0  <DIR-INPUT>
Es. $0 input_81\n"
}

if [[ ( $1 == "--help" ) || $1 == "-h" || ( -z $1 )  ]]
then
        display_usage
else
        rm -f input; ln -s $1 input
fi




