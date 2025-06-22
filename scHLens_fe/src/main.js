import Vue from "vue";
import lang from 'element-ui/lib/locale/lang/en'
import locale from 'element-ui/lib/locale'
import App from "./App.vue";
import router from "./router";
import store from "./store";
import "@/assets/css/global.less";
import "@/assets/css/font.css"
import "@/utils/dialogDrag.js"
import VueScrollTo from 'vue-scrollto';


Vue.config.productionTip = false;
Vue.use(VueScrollTo);

locale.use(lang)


new Vue({
    router,
    store,
    render: (h) => h(App),
}).$mount("#app");
