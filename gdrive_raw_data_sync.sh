# script to backup a folder to the "raw_data_backup" folder on SCG
# Need to install the gdrive CLI here: https://github.com/gdrive-org/gdrive
# and configure authentication first

if [ "$#" -ne 1 ] && [ "$#" -ne 2 ]; then
    echo ""
    echo "USAGE: gdrive_backup.sh FOLDER_TO_BACKUP [optional_remote_folder_uuid]"
    echo ""
    echo "Will sync the contents of the specified local folder to the raw_data_backup"
    echo "folder in Bhatt lab google drive:"
    echo "https://drive.google.com/drive/u/1/folders/1d3Wr2HrMBAVnkYWKqbTijahK-y1Rx_hV"
    echo "You must install the and configure the gdrive CLI before using, see:"
    echo "https://github.com/gdrive-org/gdrive"
    echo ""
    echo "Specify a different remote folder base with the second argument"
    echo ""
    exit 1
fi

# base directory to upload folders into:
sync_base=1d3Wr2HrMBAVnkYWKqbTijahK-y1Rx_hV
# get info like this: gdrive info $sync_base 

# if other uuid specified
if [ "$#" -eq 2 ]; then
    sync_base=$2
    echo "Backing up to $sync_base"
    echo $(gdrive info $sync_base)
fi

infolder=$1 
infolder_base=$(basename "$infolder")

# make directory and get uuid
n=$(gdrive mkdir -p $sync_base $infolder_base)
n2=$(echo "$n" | cut -f 2 -d " ")
# sync local with remote 
gdrive sync upload $infolder $n2
