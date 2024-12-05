import './assets/css/main.css'

import { createApp } from 'vue'
import { createPinia } from 'pinia'
import VueUswds from "vue-uswds"
import { library } from "@fortawesome/fontawesome-svg-core";
import { FontAwesomeIcon } from "@fortawesome/vue-fontawesome";
import {Tabs, Tab} from 'vue3-tabs-component';

// arrow icons
import { faArrowLeft, faArrowRight } from '@fortawesome/free-solid-svg-icons';

// social icons
import { faSquareXTwitter, faFacebookSquare, faGithub, faFlickr, faYoutubeSquare, faInstagram } from "@fortawesome/free-brands-svg-icons";
library.add(faArrowLeft, faArrowRight, faSquareXTwitter, faFacebookSquare, faGithub, faFlickr, faYoutubeSquare, faInstagram);

import App from './App.vue'
import router from './router'

const app = createApp(App)

app.use(createPinia())
app.use(VueUswds)
app.use(router)
app.component("FontAwesomeIcon", FontAwesomeIcon)
app.component("tabsGroup", Tabs)
app.component("tabItem", Tab)

app.mount('#app')
