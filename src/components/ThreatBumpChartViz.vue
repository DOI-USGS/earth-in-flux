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

    // Create a reactive object to track active categories
    const activeCategories = reactive({});

    // Define custom colors
    const customColors = [
        '#0075A2', '#FAB13C', '#6D435A', '#7FD1B9', '#EF626C'
    ];

    // Not currently used
    // const toggleCategory = (category) => {
    //     activeCategories[category] = !activeCategories[category];
    // };

    const createBumpChart = async () => {
        // Load the CSV data
        const data = await d3.csv(publicPath + 'ranked_threats.csv');
        
        // Process the data
        const nestedData = d3.group(data, d => d.Habitat, d => d.ThreatName);
        const habitats = Array.from(nestedData.keys());
        const threatNames = Array.from(new Set(data.map(d => d.ThreatName)));
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
        const margin = { top: 20, right: 20, bottom: 30, left: 50 },
            width = 800 - margin.left - margin.right,
            height = 500 - margin.top - margin.bottom;

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
            .domain([1, threatNames.length])
            .range([0, height]);

        const widthScale = d3.scaleLinear()
            .domain([0, d3.max(data, d => +d.AverageThreatMetric)])
            .range([0, 20]);

        const colorScale = d3.scaleOrdinal()
            .domain(threatCategories.value)
            .range(customColors);

        const sanitizeClass = name => name ? name.replace(/[^a-zA-Z0-9]/g, '_') : 'unknown';

        const area = d3.area()
            .x(d => x(d.habitat))
            .y0(d => y(d.rank) - widthScale(d.AverageThreatMetric) / 2)
            .y1(d => y(d.rank) + widthScale(d.AverageThreatMetric) / 2)
            .curve(d3.curveMonotoneX);

        const rankData = threatNames.map(threat => {
            const ranks = habitats.map(habitat => {
            const habitatData = nestedData.get(habitat);
            const threatData = habitatData.get(threat)[0];
            return {
                habitat: habitat,
                rank: +threatData.Rank,
                category: threatData.ThreatCategory,
                AverageThreatMetric: +threatData.AverageThreatMetric,
                ThreatName: threat
            };
            });
            return {
                threat: threat,
                values: ranks,
                category: ranks[0].category
            };
        });

        svg.selectAll('.area')
            .data(rankData)
            .enter().append('path')
            .attr('class', d => `area ${sanitizeClass(d.threat)} ${sanitizeClass(d.category)}`)
            .attr('d', d => area(d.values))
            .attr('fill', d => colorScale(d.category))
            .attr("opacity", 0.8)
            .style('stroke', 'none');

        svg.selectAll('.overlay')
            .data(rankData)
            .enter().append('path')
            .attr('class', d => `overlay ${sanitizeClass(d.threat)} ${sanitizeClass(d.category)}`)
            .attr('d', d => area(d.values))
            .attr('fill', 'transparent')
            .on('mouseover', function (event, d) {
            d3.selectAll('.area').style('opacity', 0.1);
                d3.selectAll('.point').style('opacity', 0.1);
                d3.selectAll('.label').style('opacity', 0.1);
                d3.selectAll(`.${sanitizeClass(d.threat)}`).style('opacity', 1).style('font-weight', 'bold');
            })
            .on('mouseout', function () {
                d3.selectAll('.area').style('opacity', 0.8);
                d3.selectAll('.point').style('opacity', 1);
                d3.selectAll('.label').style('opacity', 1).style('font-weight', 'normal');
            });

        svg.selectAll('.point')
            .data(rankData.flatMap(d => d.values))
            .enter().append('circle')
            .attr('class', d => `point ${sanitizeClass(d.ThreatName)} ${sanitizeClass(d.category)}`)
            .attr('cx', d => x(d.habitat))
            .attr('cy', d => y(d.rank))
            .attr('r', d => widthScale(d.AverageThreatMetric) / 2)
            .attr('fill', 'white')
            .attr('stroke', d => colorScale(d.category))
            .style('stroke-width', 2);

        svg.selectAll('.label')
            .data(rankData)
            .enter().append('text')
            .attr('class', d => `label ${sanitizeClass(d.threat)} ${sanitizeClass(d.category)}`)
            .attr('x', 30)
            .attr('y', d => y(d.values[0].rank))
            .attr('dy', '.35em')
            .attr('text-anchor', 'end')
            .attr('id', d => `label-${sanitizeClass(d.threat)}`)
            .text(d => d.threat)
            .style('font-size', '12px')
            .style('fill', d => colorScale(d.category))
            .on('mouseover', function (event, d) {
                d3.selectAll('.area').style('opacity', 0.1);
                d3.selectAll('.point').style('opacity', 0.1);
                d3.selectAll('.label').style('opacity', 0.1).style('font-weight', 'normal');
                d3.selectAll(`.${sanitizeClass(d.threat)}`).style('opacity', 1).style('font-weight', 'bold');
            })
            .on('mouseout', function () {
                d3.selectAll('.area').style('opacity', 0.8);
                d3.selectAll('.point').style('opacity', 1);
                d3.selectAll('.label').style('opacity', 1).style('font-weight', 'normal');
            });

        svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${height+5})`)
            .call(d3.axisBottom(x))
            .attr("stroke-width", 0)
            .attr("font-size", 16)
            .attr('font-weight', 700);

        const updateCategoryStyles = () => {
            threatCategories.value.forEach(category => {
            if (activeCategories[category]) {
                d3.selectAll(`.${sanitizeClass(category)}`)
                    .style('opacity', 1)
                    .style('fill', colorScale(category))
                    .style('stroke', colorScale(category));
            } else {
                d3.selectAll(`.${sanitizeClass(category)}`)
                    .style('opacity', 0.3)
                    .style('fill', 'grey20')
                    .style('stroke', 'grey20');
            }
            });
        };
};

    onMounted(() => {
        createBumpChart();
    });
</script>

<style scoped>
    .line {
        fill: none;
        stroke-width: 2px;
    }
    .point {
        fill: white;
        stroke-width: 2px;
    }
    .area {
        fill-opacity: 0.7;
    }
    .label {
        cursor: pointer;
    }
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
        /* background-color: #f8f9fa;  */
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