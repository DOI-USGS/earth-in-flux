<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
            <p v-html="text.paragraph1"/>
        </template>
        <template #figures>
            <div
                id="threat-svg-container"
                class="single maxWidth"
            >
                <ThreatsSVG 
                    class="threat-svg"
                    v-if="currentCategory == 'none'"
                />
                <HabitatThreatsSVG 
                    v-if="currentCategory == 'habitat'"
                />
                <PollutionThreatsSVG 
                    v-if="currentCategory == 'pollution'"
                />
                <ClimateThreatsSVG 
                    v-if="currentCategory == 'climate-and-weather'"
                />
                <InvasivesThreatsSVG 
                    v-if="currentCategory == 'invasive-species'"
                />
                <FishingThreatsSVG 
                    v-if="currentCategory == 'fishing-pressure'"
                />
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <div
                id="primary-icon-container"
            >
                <button
                    id="reset-button"
                    @click="currentCategory = 'none'"
                >
                    reset
                </button>
                <div
                    v-for="icon in text.iconData" :key="icon.iconID"
                >
                    <button
                        class="icon-button primary"
                        :class="{ 'active': currentCategory == icon.iconID }"
                        @click="currentCategory = icon.iconID"
                    >
                        <img class="icon-image" :src="getImageURL(icon.icon)" alt="">
                    </button>
                </div>
            </div>
            <div
                class="definition-container"
                v-for="icon in text.iconData" :key="icon.iconID" v-show="currentCategory == icon.iconID"
            >
                <p>
                    {{ icon.iconText }}
                </p>
            </div>
            <div
                v-for="icon, index in text.iconData" :key="index" v-show="currentCategory == icon.iconID"
            >
                <div
                    class="sub-icon-container"
                    v-for="subicon, index in icon.subThreatData"
                    :key="index"
                >
                    <button
                        class="icon-button"
                    >
                        <img class="icon-image" :src="getImageURL(subicon.subThreatIcon)" alt="">
                    </button>
                    <p> <span class="emph"> {{ subicon.subThreat }}</span>: {{ subicon.subThreatText }}</p>
                </div>

            </div>
        </template>
    </VizSection>
</template>

<script setup>
    import { onMounted, ref } from "vue";
    import VizSection from '@/components/VizSection.vue';
    import ThreatsSVG from "@/assets/svgs/InlandFisheriesThreats.svg";
    import HabitatThreatsSVG from "@/assets/svgs/InlandFisheriesThreats_habitat.svg";
    import PollutionThreatsSVG from "@/assets/svgs/InlandFisheriesThreats_pollution.svg";
    import ClimateThreatsSVG from "@/assets/svgs/InlandFisheriesThreats_climate-and-weather.svg";
    import InvasivesThreatsSVG from "@/assets/svgs/InlandFisheriesThreats_invasive-species.svg";
    import FishingThreatsSVG from "@/assets/svgs/InlandFisheriesThreats_fishing-pressure.svg";

    // define props
    defineProps({
        text: { type: Object }
    })

    // Global variables 
    const currentCategory = ref('none')

    onMounted(async () => {

    });

    function getImageURL(filename) {
        return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
    }
</script>

<style lang="scss" scoped>
#primary-icon-container {
    display: flex;
    gap: 30px;
    align-items: center;
    justify-content: center;
}
.sub-icon-container {
    display: flex;
    gap: 15px;
    align-items: center;
    margin-bottom: 15px;
}
#reset-button {
    background-color: transparent;
    border: 0.5px solid var(--background-color);
    box-shadow: 0px 0px 8px rgba(39,44,49,.15), 1px 4px 4px rgba(39,44,49,.04);
    padding: 0.75rem;
    border-radius: 3px;
    opacity: 0.8;
}
#reset-button:hover {
    opacity: 1;
    box-shadow: 0 0px 6px rgba(0, 0, 0, 0.2);
    transition: all .3s ease; 
}
.icon-button {
    background-color: transparent;
    border: 0.25px solid rgb(223, 223, 223);
    border-radius: 50%;
    height: 85px;
    width: 85px;
}
.primary {
    opacity: 0.5;
    border: 0.5px solid var(--background-color);
    box-shadow: 0px 0px 8px rgba(39,44,49,.15), 1px 4px 4px rgba(39,44,49,.04);
}
.active {
    opacity: 1;
    transform: scale(1.1);
}
.primary:hover {
    transform: scale(1.1);
    opacity: 1;
    transition: all .3s ease; 
}
.icon-image {
    padding: 10px;
    max-width: 75px;
    max-height: 75px;
    height: auto;
    width: auto;
    @media only screen and (max-width: 600px) {
        max-width: 40px;
        max-height: 40px;
    }
}
.definition-container {
    margin: 3rem 0 3rem 0;
}
</style>