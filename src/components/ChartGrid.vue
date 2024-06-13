<template>
    <div id="chartGrid" class="padded">
        <a v-for="(item, index) in chartContent" :key="index" :href="item.route">
            <ChartCard 
				:src="item.img_src"
				:alt="item.alt"
			/>
        </a>
    </div>
</template>

<script setup>
	import ChartCard from '@/components/ChartCard.vue'
	import ChartGrid from '@/assets/content/ChartGrid.js';
	
	const chartContent = ChartGrid.chartGridItems;

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