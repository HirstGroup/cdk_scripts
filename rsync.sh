# to copy from cluster
rsync -rvP --size-only --include='*/' --include='*strip.nc' --exclude='_*' --exclude='*.nc' --exclude='reference.frc'

# to copy to run MDs to cluster
rsync -rvP --size-only --include='*/' --include='*.parm7' --include='*.rst7' --exclude='*'