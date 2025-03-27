# scHLens - A Web Server for Hierarchically and Interactively Exploring Single Cell RNA-seq Data

This is the repository for scHLens. scHLens is an interactive analysis platform for hierarchically and interactively exploring cellular heterogeneity of scRNA-seq dataset, with various visualization views. scHLens focused on the cell type annotation of scRNA-seq. 

![system](./system.png)

## Quick start

This web-server is based on `python 3.7.12`, `R 4.2.2`, `Node.js 16.19.0`. You can launch this web-server using the following steps:

(This web-server can run on both Windows and Linux platforms. we will take Linux platform as an example.)

**Step1. Python prerequisites **

It is recommended to use Anaconda3 for managing the environment:

```
conda create --name scHLens python=3.7.13
pip install -r ./scHLens_be/requirements.txt
```

**Step2. R prerequisites**

Install `R 4.2.2`, and then run the following commands to install packages and set R_HOME:

```
Rscript ./scHLens_be/package.R
ENV R_HOME=/usr/lib/R
```

**Step3. Node.js prerequisites**

Install `Node.js 16.19.0`，and then run the following commands:

```
cd ./scHLens_fe
npm i
```

**Step4. Start front end**

```
cd ./scHLens_fe
npm start
```

**Step4. Start backend end**

```
cd ./scHLens_be
conda activate scHLens
python app.py
```

Then, you can access scHLens via http://localhost:8081



## Additional Samples





## Docker version

We also provide a more docker version 
