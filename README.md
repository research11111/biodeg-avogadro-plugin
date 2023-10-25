# Install
```bash
dir="$HOME/Library/Application Support/OpenChemistry/Avogadro/commands/biodeg" && \
 mkdir "$dir" && \
 for f in $(find . -type f -mindepth 1 -maxdepth 1) ; do \
   cp "$f" "$dir" ; \
 done
```
# Develop
```bash
dir="$HOME/Library/Application Support/OpenChemistry/Avogadro/commands/biodeg"
ln -s "$(pwd)" "$dir"
```
