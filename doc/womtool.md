Download `womtool-86.jar` from https://github.com/broadinstitute/cromwell/releases

### Set up

In the `.bashrc` file, make `womtool` executable from any directory

```bash
alias womtool='java -jar /path/to/womtool-86.jar'
```

### Validate

```bash
womtool validate workflow.wdl
```

### Visualize

Install `dot` from graphviz

```bash
sudo apt install graphviz
```

Visualize the workflow as a graph in PNG format

```bash
womtool graph workflow.wdl | dot -Tpng > workflow.png
```
