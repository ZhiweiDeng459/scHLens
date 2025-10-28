import Vue from "vue";
import VueRouter from "vue-router";

Vue.use(VueRouter);

import System from "@/System.vue"
import Home from "@/Home.vue"
import Contact from "@/Contact.vue"
import Tutorial from "@/Tutorial.vue"

import TIntroduction from "@/tutorials/TIntroduction.vue"
import TenterJob from "@/tutorials/TenterJob.vue"
import TselectDataset from "@/tutorials/TselectDataset.vue"
import TStartaAnalysisPipeline from "@/tutorials/TStartaAnalysisPipeline.vue"
import TAnalysiswithViews from "@/tutorials/TAnalysiswithViews.vue"
import THierarchicalExploration from "@/tutorials/THierarchicalExploration.vue"
import THelp from "@/tutorials/THelp.vue"
import TCase1 from "@/tutorials/TCase1.vue"
import TCase2 from "@/tutorials/TCase2.vue"
import TCase3 from "@/tutorials/TCase3.vue"

const routes = [
    {
        path:'/',
        component:Home,
    },
    {
        path:'/home',
        component:Home,
    },
    {
        path:'/system',
        component:System,
    },
    {
        path:'/tutorial',
        component:Tutorial, 
        children:[
            {
                path:'intro',
                component:TIntroduction,
            },
            {
                path:'enterjob',
                component:TenterJob,
                children:[
                    {
                        path:'1',
                        component:TenterJob,    
                    },
                    {
                        path:'2',
                        component:TenterJob,    
                    },
                    {
                        path:'3',
                        component:TenterJob,    
                    },
                    {
                        path:'4',
                        component:TenterJob,    
                    },
                    {
                        path:'5',
                        component:TenterJob,    
                    },
                    {
                        path:'6',
                        component:TenterJob,    
                    },
                ]
            },
            {
                path:'selectDataset',
                component:TselectDataset,
                children:[
                    {
                        path:'1',
                        component:TselectDataset,    
                    },
                    {
                        path:'2',
                        component:TselectDataset,    
                    },
                ]
            },
            {
                path:'startaAnalysisPipeline',
                component:TStartaAnalysisPipeline,
                children:[
                    {
                        path:'1',
                        component:TStartaAnalysisPipeline,
                    },
                    {
                        path:'2',
                        component:TStartaAnalysisPipeline,
                        children:[
                            {
                                path:'t_tSNE_perplexity',
                                component:TStartaAnalysisPipeline,
                            },
                            {
                                path:'t_UMAP_minDist',
                                component:TStartaAnalysisPipeline,    
                            },
                            {
                                path:'t_UMAP_nNeighbors',
                                component:TStartaAnalysisPipeline,    
                            },
                            {
                                path:'t_Leiden_nNeighbors',
                                component:TStartaAnalysisPipeline,    
                            },
                            {
                                path:'t_Leiden_resolution',
                                component:TStartaAnalysisPipeline,    
                            },
                        ]
                    },
                ]
            },
            {
                path:'analysiswithViews',
                component:TAnalysiswithViews,
                children:[
                    {
                        path:'1',
                        component:TAnalysiswithViews,
                        children:[
                            {
                                path:'1',
                                component:TAnalysiswithViews,
                            },
                            {
                                path:'2',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'3',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'4',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'5',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'6',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'7',
                                component:TAnalysiswithViews,    
                            },
                        ]
        
                    },
                    {
                        path:'2',
                        component:TAnalysiswithViews, 
                        children:[
                            {
                                path:'1',
                                component:TAnalysiswithViews,
                            },
                            {
                                path:'2',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'3',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'4',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'5',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'6',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'7',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'8',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'9',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'10',
                                component:TAnalysiswithViews,    
                            },
                            {
                                path:'11',
                                component:TAnalysiswithViews,    
                            },
                        ]
  
                    },
                ]

            },
            {
                path:'hierarchicalExploration',
                component:THierarchicalExploration,
                children:[
                    {
                        path:'1',
                        component:THierarchicalExploration,
                    },
                    {
                        path:'2',
                        component:THierarchicalExploration,    
                    },
                    {
                        path:'3',
                        component:THierarchicalExploration,    
                    },
                    {
                        path:'4',
                        component:THierarchicalExploration,    
                    },
                ]
            },
            {
                path:'help',
                component:THelp,
            },       
            {
                path:'case1',
                component:TCase1,
            },

            {
                path:'case2',
                component:TCase2,
            },            {
                path:'case3',
                component:TCase3,
            },
        ],
    },
    {
        path:'/contact',
        component:Contact,
    },
];

//忽略重复使用一个url进行导航时报错
const originPush = VueRouter.prototype.push
VueRouter.prototype.push = function push(path){
    return originPush.call(this,path).catch(err=>err)
}

const router = new VueRouter({
    routes,
});

export default router;
