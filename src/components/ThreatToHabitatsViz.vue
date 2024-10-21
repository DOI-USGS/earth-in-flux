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
                <p v-html="text.paragraph1" />
        </template>
        <template #figures>
            <div id="threat-container" ref="chart"></div>
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

    // define props
    defineProps({
        text: { type: Object }
    })

    // Global variables 
    const publicPath = import.meta.env.BASE_URL;
    const mobileView = isMobile;
    const data = ref();
    const dataFile = 'findex_total_weighted_threats.csv'
    const chart = ref(null);
    let chartDimensions;
    let chartBounds;
    let nodeGroup;
    let linkGroup;
    let textGroup;
    const nodeWidth = 4;
    const labelBuffer = 10;

    // Colors for threat categories, Needs to be updated with CSS for text legend
    const categoryColors = {
        'Climate and weather': '#c29fcd',
        'Exploitation':  '#d38884',
        'Habitat':  '#dab589', 
        'Invasive species':  '#729C9D',
        'Pollution': '#899bb7'
    }; 

    onMounted(async () => {
        try {
            await loadDatasets();
            
            if (data.value.length > 0) {
                initLineChart({
                    width: chart.value.offsetWidth,
                    height: window.innerHeight * 0.8,
                    margin: 10,
                    marginLeft: mobileView ? 80: 150,
                    marginRight: mobileView ? 125: 250,
                    marginTop: 30,
                    containerId: 'threat-container'
                });
                createLineChart({
                    dataset: data.value,
                    containerId: 'threat-container'
                });
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

    function initLineChart({
        width,
        height,
        margin,
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin, // right margin, in pixels
        containerId
    }) {
        // set up global chart dimensions, including bounded dimensions
        chartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            },
            boundedWidth: width - marginLeft - marginRight,
            boundedHeight: height - marginTop - marginBottom
        }


        // draw svg canvas for LineChart
        const svg = d3.select('#' + containerId)
            .append('svg')
            .attr('class', 'LineChartSVG')
            .attr('viewBox', `0 0 ${chartDimensions.width} ${chartDimensions.height}`)
            .style('width', "100%")
            .style('height', "100%");

        // add group for bar chart bounds, translating by chart margins
        chartBounds = svg.append('g')
            .attr('id', 'wrapper')
            .style("transform", `translate(${
                chartDimensions.margin.left
            }px, ${
                chartDimensions.margin.top
            }px)`)


    };

    function createLineChart({
        dataset
    }) {

        
    };

    
</script>

<style lang="scss">
    #threat-container {
        max-width: 1000px;
        margin: 5rem auto 0 auto;
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
</style>