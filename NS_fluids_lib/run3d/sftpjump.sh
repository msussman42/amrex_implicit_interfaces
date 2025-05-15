# use get or put
# alternative: rsync -rv <source> <dest>
# alternative: rsync -rv -e 'ssh -J sussman@henri.math.fsu.edu:22'  sussman@compute1.math.fsu.edu:/scratch1/sussman/amrex_implicit_interfaces/NS_fluids_lib/run3d/tank_sloshing/temp.tar.gz .
sftp -oProxyJump=sussman@henri.math.fsu.edu sussman@bromwich.math.fsu.edu
