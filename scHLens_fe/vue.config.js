module.exports = {
    publicPath:'./',
    devServer: {
        host: "0.0.0.0",
        port: "8081",
        filenameHashing: true,
        // https: true,
        proxy: {
            "/api": { //普通api
                target: "http://127.0.0.1:5003/scHLens",
                changeOrigin: true,
                pathRewrite: {
                    "": "",
                },
            },
            "/socket.io": { //socket请求
                target: "http://127.0.0.1:5003",
                changeOrigin: true,
                pathRewrite: {
                    "": "",
                },
            },
        },
    },
    configureWebpack: {
        devtool: 'source-map'
    }
    
};




