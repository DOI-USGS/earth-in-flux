<template>
    <div id="chartGrid" class="padded">
		<ChartCard @click.enter="showSubPage(item.project, item.vizRoute)" v-for="(item, index) in filteredChartContent" :key="index"
			:id="item.vizRoute"
			:src="getThumb(item.img_src)"
			:alt="item.alt"
		/>
    </div>
</template>

<script setup>
	import { useRouter } from 'vue-router';
	import { computed } from 'vue';

	import ChartCard from '@/components/ChartCard.vue';
	import ChartGrid from '@/assets/content/ChartGrid.js';

	const props = defineProps({
        view: {
            type: String,
            default: ``
        },
    })
	
	// global variables
	const router = useRouter()
	const chartContent = ChartGrid.chartGridItems;
	// set up filtered chart data as computed property
	const filteredChartContent = computed(() => {
		return props.view == 'all' ? chartContent : chartContent.filter(d => d.project.replace(/\s+/g, '-').toLowerCase() === props.view)
	});

	function showSubPage(project, vizRoute) {
		const projectRoute = project.replace(/\s+/g, '-').toLowerCase();
		router.push({ name: 'SubPage', params: { projectRoute, vizRoute } })
	}

	function getThumb(pic) {
        return 'https://labs.waterdata.usgs.gov/visualizations/thumbnails/'+pic
    }

</script>

<style>

#chartGrid{
	margin: 50px auto;
	display: flex;
	justify-content: center;
	gap: 50px;
	flex-wrap: wrap;
	gap: 50px;
	width: 100%;
	max-width: 60%;
}
</style>