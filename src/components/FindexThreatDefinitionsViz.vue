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
            >
                <component 
                    :is="getCurrentSVG()" 
                    class="threat-svg"
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
                    v-for="icon in text.iconData"
                    :key="icon.iconID"
                    class="button-wrapper"
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
                    <div
                        class="button-wrapper-vertical"
                    >
                        <button
                            class="icon-button sub"
                            :class="icon.iconID"
                        >
                            <img class="icon-image" :src="getImageURL(subicon.subThreatIcon)" alt="">
                        </button>
                    </div>
                    <p> <span class="emph"> {{ subicon.subThreat }}</span>: {{ subicon.subThreatText }}</p>
                </div>

            </div>
        </template>
    </VizSection>
</template>

<script setup>
    import { ref } from "vue";
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

    const getCurrentSVG = () => {
        const components = {
            none: ThreatsSVG,
            habitat: HabitatThreatsSVG,
            pollution: PollutionThreatsSVG,
            'climate-and-weather': ClimateThreatsSVG,
            'invasive-species': InvasivesThreatsSVG,
            'fishing-pressure': FishingThreatsSVG
        };
        return components[currentCategory.value] || ThreatsSVG;
    }

    function getImageURL(filename) {
        return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
    }
</script>

<style lang="scss" scoped>
#threat-svg-container {
    display: flex;
    justify-content: center;
}
.threat-svg {
    width: 100%;
    height: 100%;
    max-width: 700px;
}
#primary-icon-container {
    display: flex;
    align-items: center;
    justify-content: center;
    max-width: 100%;
}
.sub-icon-container {
    display: flex;
    gap: 15px;
    align-items: center;
    margin-bottom: 15px;
}
.sub-icon-container p {
    padding: 0;
}
#reset-button {
    background-color: transparent;
    border: 0.5px solid var(--background-color);
    box-shadow: 0px 0px 8px rgba(39,44,49,.15), 1px 4px 4px rgba(39,44,49,.04);
    padding: 0.75rem;
    border-radius: 3px;
    opacity: 0.8;
    margin: 3px;
}
#reset-button:hover {
    opacity: 1;
    box-shadow: 0 0px 6px rgba(0, 0, 0, 0.2);
    transition: all .3s ease; 
}
.button-wrapper {
    display: flex;
    flex: 1 1 50%;
    max-width: 115px;
}
.button-wrapper-vertical {
    display: flex;
    flex: 1 1 50%;
    max-width: 120px;
    @media screen and (max-width: 600px) {
        max-width: 60px;
    }
}
.icon-button {
    flex: 1;
    background-color: transparent;
    border: 0.25px solid rgb(223, 223, 223);
    border-radius: 50%;
    padding: 0px;
    margin: 15px;
    @media screen and (max-width: 600px) {
        margin: 5px;
    }
}
.primary {
    min-width: 35px;
    opacity: 0.5;
    border: 0.5px solid var(--background-color);
    box-shadow: 0px 0px 8px rgba(39,44,49,.15), 1px 4px 4px rgba(39,44,49,.04);
    @media screen and (max-width: 600px) {
        margin: 7px;
    }
}
.primary:hover {
    transform: scale(1.1);
    opacity: 1;
    transition: all .3s ease; 
}
.sub {
    min-width: 55px;
}
.icon-button:after {
    content: '';
    display: inline-block;
    vertical-align: top;
    padding-top: 100%;
}
.active {
    opacity: 1;
    transform: scale(1.1);
}
.habitat {
    border-color: var(--color-habitat-faded);
}
.pollution {
    border-color: var(--color-pollution-faded);
}
.climate-and-weather {
    border-color: var(--color-climate-and-weather-faded);
}
.icon-image {
    transform: translate(0, 32%); /* negate vertical-align top on parent*/
    max-height: 60%;
    max-width: 60%;
}
.definition-container {
    margin: 3rem 0 3rem 0;
}
</style>