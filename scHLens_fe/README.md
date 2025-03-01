# sc-projection

## Project setup
```
yarn install
```

### Compiles and hot-reloads for development
```
yarn start
```

### Compiles and minifies for production
```
yarn build
```

### Lints and fixes files
```
yarn lint
```

### Customize configuration
See [Configuration Reference](https://cli.vuejs.org/config/).

### QA

1. 为什么D3.select会查到Vue文件的名字(or 组件name)？例如d3.select('.global-scatter')会匹配GlobalScatter.vue的根元素？
2. 为什么second-row不加overflow的话内容会撑开center这个flexbox?