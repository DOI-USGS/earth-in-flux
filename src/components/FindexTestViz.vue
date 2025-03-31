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
                v-model="selected"
                :options="[
                    { label: 'Climate and weather', value: 'Climate_and_weather_map', color: '#2a9d8f' },
                    { label: 'Fishing pressure', value: 'Fishing_pressure_map', color: '#264653' },
                    { label: 'Habitat', value: 'Habitat_map', color: '#899bb7' },
                    { label: 'Invasive species', value: 'Invasive_species_map', color: '#c29fcd' },
                    { label: 'Pollution', value: 'Pollution_map', color: '#dab589' }
                ]"
            />
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div id="image-container" ref="chart">
                <img
                    :src="`@/assets/images/${selected}.png`"
                    :alt="selected"
                    class="preview-img"
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

    // define props
    defineProps({
        text: { type: Object }
    })

    // Global variables 
    const publicPath = import.meta.env.BASE_URL;
    const mobileView = isMobile;
    const data = ref();
    //const dataFile = 'findex_total_weighted_threats.csv'
    const chart = ref(null);
    let chartDimensions;
    let chartBounds;

    const selected = ref('mountains')


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

    function initSankey({
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


        // draw svg canvas for sankey
        const svg = d3.select('#' + containerId)
            .append('svg')
            .attr('class', 'sankeySVG')
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

        // Add group to chart bounds to hold all sankey path groups
        nodeGroup = chartBounds.append('g')
            .attr('id', 'node_group')
        
        linkGroup = chartBounds.append('g')
            .attr('id', 'link_group')

        textGroup = chartBounds.append('g')
            .attr('id', 'text_group')

        // add titles
        const subtitle = mobileView ? "High to low" : "Ranked high to low"
        const leftTitle = svg.append("text")
            .attr("class", "axis-title")
            .attr("x", chartDimensions.margin.left - labelBuffer + nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", "0em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.left)
            .style("text-anchor", "end")
            .text("Threat Categories")
            .call(d => mobileView ? wrap(d, {shift: false}) : d)

        const leftTitleLength = leftTitle.node().getComputedTextLength()

        svg.append("text")
            .attr("class", "axis-text axis-value axis-notation")
            .attr("x", chartDimensions.margin.left - labelBuffer + nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", leftTitleLength > chartDimensions.margin.left ? "2.8em" : "1.5em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.left)
            .attr("text-anchor", "end")
            .text(subtitle)
            .call(d => wrap(d, {shift: false}))

        const rightTitle = svg.append("text")
            .attr("class", "axis-title")
            .attr("x", chartDimensions.width - chartDimensions.margin.right + labelBuffer - nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", "0em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.right)
            .style("text-anchor", "start")
            .text("Threats")

        const rightTitleLength = rightTitle.node().getComputedTextLength()

        svg.append("text")
            .attr("class", "axis-text axis-value axis-notation")
            .attr("x", chartDimensions.width - chartDimensions.margin.right + labelBuffer - nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", rightTitleLength > chartDimensions.margin.right ? "2.8em" : "1.5em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.right)
            .attr("text-anchor", "start")
            .text(subtitle)
            .call(d => wrap(d, {shift: false}))
    };


</script>

<style lang="scss">
    #threat-container {
        max-width: min(1000px, 60vw);
        margin: 5rem auto 0 auto;
        @media screen and (max-width: 600px) {
            max-width: 100%;
        }
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