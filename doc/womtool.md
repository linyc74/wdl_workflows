Download `womtool-86.jar` from https://github.com/broadinstitute/cromwell/releases

### Validate

```bash
java -jar womtool-86.jar validate SomaticPipelineTumorNormalMode.wdl
```

### Visualize

Install `dot` from graphviz

```bash
sudo apt install graphviz
```

Visualize the workflow as a graph in PNG format

```bash
java -jar womtool-86.jar graph SomaticPipelineTumorNormalMode.wdl | dot -Tpng > SomaticPipelineTumorNormalMode.png
```
