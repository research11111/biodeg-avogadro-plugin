# Install
## MacOS
```bash
dir="$HOME/Library/Application Support/OpenChemistry/Avogadro/commands/biodeg" && \
 mkdir "$dir" && \
 for f in $(find . -type f -mindepth 1 -maxdepth 1) ; do \
   cp "$f" "$dir" ; \
 done
```
## Ubuntu
```bash
dir="$HOME/.local/share/OpenChemistry/Avogadro/commands/biodeg" && \
 mkdir "$dir" && \
 for f in $(find . -type f -mindepth 1 -maxdepth 1) ; do \
   cp "$f" "$dir" ; \
 done
```
# Develop
## MacOS
```bash
dir="$HOME/Library/Application Support/OpenChemistry/Avogadro/commands/biodeg"
ln -s "$(pwd)" "$dir"
```
## Ubuntu
```bash
dir="$HOME/.local/share/OpenChemistry/Avogadro/commands/biodeg"
ln -s "$(pwd)" "$dir"
```
