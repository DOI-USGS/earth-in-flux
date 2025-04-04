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
            <RadioGroup
                v-model="selectedLayer"
                :options="layers"
            />
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div id="image-container" ref="chart">
                <RegionMap
                    :selected-layer="selectedLayer"
                    :layer-paths="layerPaths"
                    :layer-x="layerX"
                    :layer-y="layerY"
                    :layer-mag="1"
                    :topo-regions="topoRegions"
                    :regions-var="Region_nam"
                    :topo-us="topoUS"
                />
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
        </template>
    </VizSection>
</template>

<script setup>
    import { onMounted, ref } from "vue";
    import { isMobile } from 'mobile-device-detect';
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';
    import RadioGroup from './RadioGroup.vue'
    import RegionMap from './RegionMap.vue'

    import ClimateMap from '@/assets/images/Climate_and_weather_map.png'
    import FishingMap from '@/assets/images/Fishing_pressure_map.png'
    import HabitatMap from '@/assets/images/Habitat_map.png'
    import InvasiveMap from '@/assets/images/Invasive_species_map.png'
    import PollutionMap from '@/assets/images/Pollution_map.png'
    import topoRegions from '@/assets/geo/Regions.json'
    import topoUS from '@/assets/geo/USoutline.json'


    const layers = [
        {
            label: 'Climate and weather',
            value: 'Climate_and_weather_map',
            color: '#2a9d8f',
            path: ClimateMap
        },
        {
            label: 'Fishing pressure',
            value: 'Fishing_pressure_map',
            color: '#264653',
            path: FishingMap
        },
        {
            label: 'Habitat',
            value: 'Habitat_map',
            color: '#899bb7',
            path: HabitatMap
        },
        {
            label: 'Invasive species',
            value: 'Invasive_species_map',
            color: '#c29fcd',
            path: InvasiveMap
        },
        {
            label: 'Pollution',
            value: 'Pollution_map',
            color: '#dab589',
            path: PollutionMap
        }
    ]
    const layerPaths = Object.fromEntries(
        layers.map(layer => [layer.value, { path: layer.path }])
    )



    // define props
    defineProps({
        text: { type: Object }
    })

    // Global variables 
    const publicPath = import.meta.env.BASE_URL;
    const data = ref([]);
    const chart = ref(null);

    const selectedLayer = ref('Climate_and_weather_map')


    onMounted(async () => {
        try {
            //await loadDatasets();
            
            if (data.value.length > 0) {
                
                
            } else {
                console.error('Error loading data');
            }
        } catch (error) {
            console.error('Error during component mounting', error);
        }
    });

    async function loadDatasets() {
        try {
            data.value = await loadData(dataFile);
            console.log('data in');
        } catch (error) {
            console.error('Error loading datasets', error);
        }
    };

    async function loadData(fileName) {
        try {
            const data = await d3.csv(publicPath + fileName, d => {
                return d;
            });
            return data;
        } catch (error) {
            console.error(`Error loading data from ${fileName}`, error);
            return [];
        }
    };

</script>

<style lang="scss">
    #image-container {
        max-width: min(1400px, 60vw);
        margin: 5rem auto 0 auto;
        @media screen and (max-width: 600px) {
            max-width: 100%;
        }
    }
    .preview-img {
        max-width: 100%;
        height: auto;
    }
    .axis-text {
        font-size: 1.6rem;
        font-family: var(--default-font);
        user-select: none;
        @media screen and (max-width: 600px) {
            font-size: 1.4rem;
        }
    }
    .axis-title {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 900;
        fill: var(--color-text);
        user-select: none;
        @media screen and (max-width: 600px) {
            font-size: 1.6rem;
        }
    }
    .axis-notation {
        font-style: italic;
    }
    
</style>