/etc
/var/adm
/var/log

ls -ultar
ls -cltar

sudo find . -mtime -0 -exec touch {} \;

find ./ -print | while read filename; do
    # do whatever you want with the file
    touch -d "$(date -R -r "$filename") - 2 hours" "$filename"
done
