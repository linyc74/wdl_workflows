Download `womtool-86.jar` from https://github.com/broadinstitute/cromwell/releases

### Validate

```bash
java -jar womtool-86.jar validate workflow.wdl
```

### Visualize

Install `dot` from graphviz

```bash
sudo apt install graphviz
```

Visualize the workflow as a graph in PNG format

```bash
java -jar womtool-86.jar graph workflow.wdl | dot -Tpng > workflow.png
```
