<template>
    <div id="chartGrid" class="padded">
		<ChartCard @click.enter="showSubPage(item.project, item.route)" v-for="(item, index) in chartContent" :key="index"
			:src="item.img_src"
			:alt="item.alt"
		/>
    </div>
</template>

<script setup>
	import { useRouter } from 'vue-router'

	import ChartCard from '@/components/ChartCard.vue'
	import ChartGrid from '@/assets/content/ChartGrid.js';
	
	// global variables
	const router = useRouter()
	const chartContent = ChartGrid.chartGridItems;

	function showSubPage(project, route) {
		const projectRoute = project.replace(/\s+/g, '-').toLowerCase();
		router.push(`${projectRoute}/${route}`)
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