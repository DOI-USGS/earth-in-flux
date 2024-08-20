<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
                {{ text.heading }}
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="chart-container" ref="chart">
                <div class="toggles">
                    <span v-for="category in threatCategories" :key="category">
                        <input type="checkbox" :id="category" v-model="activeCategories[category]" />
                        <label :for="category">{{ category }}</label>
                    </span>
                </div>
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
    import { onMounted, ref, reactive, watch } from 'vue';
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // global variables
    const publicPath = import.meta.env.BASE_URL;
    const chart = ref(null);
    const threatCategories = ref([]);
    const defaultOpacity = 0.8;
    const dimOpacity = 0.15;

    // Create a reactive object to track active categories
    const activeCategories = reactive({});

    // Define custom colors
    const customColors = [
        '#0075A2', '#FAB13C', '#6D435A', '#7FD1B9', '#EF626C'
    ];

    // Behavior on mounted - functions called here
    onMounted(() => {
        createBumpChart();
    });

    const createBumpChart = async () => {
        // Load the CSV data
        const data = await d3.csv(publicPath + 'findex_ranked_threats.csv');
        
        // Process the data
        const nestedData = d3.group(data, d => d.Habitat, d => d.ThreatCode);
        const habitats = Array.from(nestedData.keys());
        const threatCodes = Array.from(new Set(data.map(d => d.ThreatCode)));
        threatCategories.value = Array.from(new Set(data.map(d => d.ThreatCategory)));
  
        // Initialize the activeCategories object
        threatCategories.value.forEach(category => {
            activeCategories[category] = true;
        });

        // Watch for changes in active categories and update the chart
        watch(activeCategories, () => {
            updateCategoryStyles();
        }, { deep: true });

        // Set up the SVG canvas dimensions
        const margin = { top: 20, right: 20, bottom: 450, left: 200 },
            width = 800 - margin.left - margin.right,
            height = 1000 - margin.top - margin.bottom;

        const svg = d3.select(chart.value)
            .append('svg')
            .attr('width', width + margin.left + margin.right)
            .attr('height', height + margin.top + margin.bottom)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // Set up scales
        const x = d3.scalePoint()
            .domain(habitats)
            .range([0, width])
            .padding(0.5);

        const y = d3.scaleLinear()
            .domain([1, threatCodes.length])
            .range([0, height]);

        const widthScale = d3.scaleLinear()
            .domain([0, d3.max(data, d => +d.AverageThreatMetric)])
            .range([10, height / threatCodes.length]);

        const colorScale = d3.scaleOrdinal()
            .domain(threatCategories.value)
            .range(customColors);

        const area = d3.area()
            .x(d => x(d.habitat))
            .y0(d => y(d.rank) - widthScale(d.AverageThreatMetric) / 2)
            .y1(d => y(d.rank) + widthScale(d.AverageThreatMetric) / 2)
            .curve(d3.curveMonotoneX);

        const rankData = threatCodes.map(threatCode => {
            const ranks = habitats.map(habitat => {
                const habitatData = nestedData.get(habitat);
                const threatData = habitatData.get(threatCode)[0];
                return {
                    habitat: habitat,
                    rank: +threatData.Rank,
                    category: threatData.ThreatCategory,
                    AverageThreatMetric: +threatData.AverageThreatMetric,
                    ThreatCode: threatCode,
                    threat: threatData.Threat
                };
            });
            return {
                threatCode: threatCode,
                values: ranks,
                category: ranks[0].category
            };
        });
        
        const bumpAreasPoints = svg.append("g")
            .attr("id", "areas-points")

        bumpAreasPoints
            .selectAll('.area')
                .data(rankData)
                .enter().append('path')
                .attr('class', d => `area ${sanitizeClass(d.threatCode)} ${sanitizeClass(d.category)}`)
                .attr('d', d => area(d.values))
                .attr('fill', d => colorScale(d.category))
                .style("opacity", defaultOpacity)
                .style('stroke', 'none');

        bumpAreasPoints
            .selectAll('.point')
                .data(rankData.flatMap(d => d.values))
                .enter().append('circle')
                .attr('class', d => `point ${sanitizeClass(d.ThreatCode)} ${sanitizeClass(d.habitat)} ${sanitizeClass(d.category)}`)
                .attr('cx', d => x(d.habitat))
                .attr('cy', d => y(d.rank))
                .attr('r', d => widthScale(d.AverageThreatMetric) / 2)
                .attr('fill', 'white')
                .attr('stroke', d => colorScale(d.category))
                .style('stroke-width', 2)
                .on('mouseover', function (event, d) {
                    mouseoverThreat(d.ThreatCode, d.category)
                })
                .on('mouseout', mouseoutThreat);

        svg.append("g")
            .attr("id", "overlays")
            .selectAll('.overlay')
                .data(rankData)
                .enter().append('path')
                .attr('class', d => `overlay ${sanitizeClass(d.threatCode)} ${sanitizeClass(d.category)}`)
                .attr('d', d => area(d.values))
                .attr('fill', 'transparent')
                .on('mouseover', function (event, d) {
                    mouseoverThreat(d.threatCode, d.category)
                })
                .on('mouseout', mouseoutThreat);
        
        svg.append("g")
            .attr("id", "labels")
            .selectAll('.label')
                .data(rankData)
                .enter().append('text')
                .attr('class', d => `label ${sanitizeClass(d.threatCode)} ${sanitizeClass(d.category)}`)
                .attr('x', 10)
                .attr('y', d => y(d.values[0].rank))
                .attr('dy', '.35em')
                .attr('text-anchor', 'end')
                .attr('id', d => `label-${sanitizeClass(d.threatCode)}`)
                .text(d => d.values[0].threat)
                .style('font-size', '12px')
                .style('fill', d => colorScale(d.category))
                .on("mouseover", function (event, d) {
                    mouseoverThreat(d.threatCode, d.category)
                })
                .on('mouseout', mouseoutThreat);

        const xAxis = svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${height + 15})`)
            .call(d3.axisBottom(x))
            .attr("stroke-width", 0)
            .attr("font-size", 16)
            .attr('font-weight', 700)

        xAxis.selectAll('text')
            .attr('text-anchor', 'end')
            .attr('dy', '-0.5rem')
            .attr('transform', 'rotate(270)');

        const updateCategoryStyles = () => {
            threatCategories.value.forEach(category => {
                if (activeCategories[category]) {
                    d3.selectAll(`.overlay.${sanitizeClass(category)}`)
                        .raise()
                    d3.selectAll(`.area.${sanitizeClass(category)}`)
                        .style('opacity', defaultOpacity)
                        .style('fill', colorScale(category));
                    d3.selectAll(`.point.${sanitizeClass(category)}`)
                        .style('opacity', 1)
                        .style('stroke', colorScale(category));
                    d3.selectAll(`.label.${sanitizeClass(category)}`)
                        .style('opacity', 1)
                        .style('fill', colorScale(category));
                } else {
                    d3.selectAll(`.overlay.${sanitizeClass(category)}`)
                        .lower()
                    d3.selectAll(`.area.${sanitizeClass(category)}`)
                        .style('opacity', dimOpacity)
                        .style('fill', '#949494');
                    d3.selectAll(`.point.${sanitizeClass(category)}`)
                        .style('opacity', dimOpacity)
                        .style('stroke', '#949494');
                    d3.selectAll(`.label.${sanitizeClass(category)}`)
                        .style('opacity', 1)
                        .style('fill', '#6E6E6E');
                }
            });
        };
    };

    

    function sanitizeClass(name) {
        return name ? name.replace(/[^a-zA-Z0-9]/g, '_') : 'unknown';
    } 

    function mouseoverThreat(threatCode, category) {
        if (activeCategories[category]) {
            d3.selectAll('.area')
                .transition(getUpdateTransition())
                .style('opacity', dimOpacity);
            d3.selectAll('.point')
                .transition(getUpdateTransition())
                .style('stroke-opacity', dimOpacity);
            d3.selectAll('.label')
                .transition(getUpdateTransition())
                .style('opacity', dimOpacity);
            d3.selectAll(`.${sanitizeClass(threatCode)}`)
                .transition(getUpdateTransition())
                .style('opacity', 1)
                .style('stroke-opacity', 1)
                .style('font-weight', 'bold');
            d3.selectAll(`.area.${sanitizeClass(threatCode)}`)
                .raise();
            d3.selectAll(`.point.${sanitizeClass(threatCode)}`)
                .raise();
        }
    }

    function mouseoutThreat() {
        threatCategories.value.forEach(category => {
            if (activeCategories[category]) {
                d3.selectAll(`.area.${sanitizeClass(category)}`)
                    .transition(getUpdateTransition())
                    .style('opacity', defaultOpacity);
                d3.selectAll(`.point.${sanitizeClass(category)}`)
                    .transition(getUpdateTransition())
                    .style('stroke-opacity', 1);
                d3.selectAll(`.label.${sanitizeClass(category)}`)
                    .transition(getUpdateTransition())
                    .style('opacity', 1)
                    .style('font-weight', 'normal');
            }
        });
        d3.selectAll('.point')
            .raise();
    }

    function getUpdateTransition() {
      return d3.transition()
        .duration(200)
        .ease(d3.easeCubicInOut)
    }

</script>

<style scoped>
    .chart-container {
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        width: 100%;
        max-width: 900px;
        height: auto;
        margin: auto;
        margin-top: 50px;
    }
    .toggles {
        display: flex;
        justify-content: center;
        align-items: center;
        margin-bottom: 20px;
    }
    .toggles span {
        margin-right: 10px;
    }
    .toggles label {
        margin-left: 5px;
    }
</style>
<style lang="scss">
    .label {
        cursor: pointer;
    }
</style>