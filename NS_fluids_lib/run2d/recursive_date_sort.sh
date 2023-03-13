find . -type f -printf "%T@ %Tc %p\n" | sort -n > out.txt
