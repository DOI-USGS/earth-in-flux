import { createRouter, createWebHistory } from 'vue-router'
import VisualizationView from '@/views/VisualizationView.vue'
import SubPage from '@/components/SubPage.vue'

function lazyLoad(view){
  return() => import(`@/views/${view}.vue`)
}

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: '/',
      name: 'VisualizationContent',
      component: VisualizationView
    },
    {
      path: '/index.html',
      name: 'Index',
      component: VisualizationView
    },
    {
      path: '/:projectRoute/:route',
      name: 'Subpage',
      component: SubPage
    },
    {
      path: "/404",
      name: "Error404",
      component: lazyLoad('Error404Page')
    },
    { 
      path: '/:pathMatch(.*)*', 
      redirect: { name: "Error404" }
    }
  ]
})

export default router
