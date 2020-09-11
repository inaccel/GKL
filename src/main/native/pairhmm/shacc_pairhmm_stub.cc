#include "shacc_pairhmm.h"
#include <debug.h>

#include <cstring>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <vector>

#include <inaccel/coral>
#include <cstdlib>
#include <iostream>
#include <time.h>

bool calculate(Batch& batch);

int read_pad_length(int length){
	return length % COLS ? COLS - length % COLS : 0;
}

int hap_pad_length(int length){
	return length < ROWS ? ROWS - length : 0;
}


void add_reads(inaccel::vector<ReadData> &read_data, Read *reads, int offset, int num_reads) {

	int index = 0;
	for (int n = offset; n < offset + num_reads; n++) {
		int read_length = reads[n].length;
		int read_length_padded = read_length + read_pad_length(read_length);

		for (int i = 0; i < read_length; i++) {
			read_data[index + i].read_num = offset + n;
			read_data[index + i].base = reads[n].bases[i];
			read_data[index + i].position = 0;
			read_data[index + i].position |= i == 0 ? 1 : 0;
			read_data[index + i].position |= i == read_length - 1 ? 2 : 0;

			float mm = 0.1f;
			float gm = 0.2f;
			read_data[index + i].mx = 0.3f;
			read_data[index + i].my = 0.4f;
			read_data[index + i].gg = 0.5f;
			read_data[index + i].mm_1m_qual = mm * 0.6f;
			read_data[index + i].mm_qual_div3 = mm * 0.7f / 3.0f;
			read_data[index + i].gm_1m_qual = gm * 0.8f;
			read_data[index + i].gm_qual_div3 = gm * 0.9f / 3.0f;
		}

		for (int i = read_length; i < read_length_padded; i++) {
			read_data[index + i].base = '-';
			read_data[index + i].position = 0;
		}
		index += read_length_padded;
	}
}

void add_haps(inaccel::vector<HapData> &hap_data, Haplotype *haps, int offset, int num_haps) {

	int index = 0;
	for (int n = offset; n < offset + num_haps; n++) {
		int hap_length = haps[n].length;
		int hap_length_padded = hap_length + hap_pad_length(hap_length);

		for (int i = 0; i < hap_length; i++) {
			hap_data[index + i].base = haps[n].bases[i];
			hap_data[index + i].y_init = 1.33e36 / hap_length;
			hap_data[index + i].hap_num = offset + n;
			hap_data[index + i].position = 0;
			hap_data[index + i].position |= i == 0 ? 1 : 0;
			hap_data[index + i].position |= i == hap_length - 1 ? 2 : 0;
		}

		for (int i = hap_length; i < hap_length_padded; i++) {
			hap_data[index + i].base = '-';
			hap_data[index + i].position = 0;
		}
		index += hap_length_padded;
	}
}

float fpga_pairhmm(testcase testcase){
	Batch m_batch;
	std::vector<Read> m_reads;
	std::vector<Haplotype> m_haps;

	m_batch.num_reads = 1;
	m_batch.num_haps = 1;

	int num_testcases = m_batch.num_reads * m_batch.num_haps;

	// get read

	Read read;
	read.bases = testcase.rs;
	read.length = testcase.rslen;
      	read.i = testcase.i;
      	read.d = testcase.d;
      	read.c = testcase.c;
      	read.q = testcase.q;
      	m_reads.push_back(read);

    	m_batch.reads = m_reads.data();

    	// get haplotype

      	Haplotype hap;
      	DBG("hap #%d len = %d", i, testcase.haplen);
      	hap.bases = testcase.hap;
      	hap.length = testcase.haplen;
      	m_haps.push_back(hap);

    	m_batch.haps = m_haps.data();

    	// allocate results
    	m_batch.results = (float*)malloc(sizeof(float) * num_testcases);

    	bool valid = calculate(m_batch);

	if (valid) {
		return m_batch.results[0];
	}

}

bool calculate(Batch& batch){
	std::vector<inaccel::vector<HapData>> haps_buffers;
	std::vector<inaccel::vector<ReadData>> reads_buffers;
	std::vector<int> reads_counts;
	std::vector<int> haps_counts;

	int read_buf = 0;

	std::vector<int> total_read_length = {0};
	std::vector<int> max_read_length = {0};

	int cnt = 0;
	// Fill inaccel vectors with Read Data
	for (int i = 0; i < batch.num_reads; i++) {
		int read_length = batch.reads[i].length;
		int read_length_padded = read_length + read_pad_length(read_length);

		if (total_read_length[read_buf] + read_length_padded >= MAX_READ_LENGTH) {
			inaccel::vector<ReadData> reads(total_read_length[read_buf]);
			add_reads(reads, batch.reads, i - cnt, cnt);

			reads_buffers.push_back(reads);

			reads_counts.push_back(cnt);
			cnt = 0;

			read_buf++;

			total_read_length.push_back(0);
			max_read_length.push_back(0);
		}

		cnt++;
		total_read_length[read_buf] += read_length_padded;
		max_read_length[read_buf] = std::max(read_length, max_read_length[read_buf]);

	}

	inaccel::vector<ReadData> reads(total_read_length[read_buf]);
	add_reads(reads, batch.reads, batch.num_reads - cnt, cnt);
	reads_buffers.push_back(reads);

	reads_counts.push_back(cnt);

	cnt = 0;
	int hap_buf = 0;
	std::vector<int> total_hap_length = {0};

	// Fill inaccel vectors with Hap Data
	for (int i = 0; i < batch.num_haps; i++) {
		int hap_length = batch.haps[i].length;
		int hap_length_padded = hap_length + hap_pad_length(hap_length);

		if (total_hap_length[hap_buf] + hap_length_padded >= MAX_HAP_LENGTH) {
			inaccel::vector<HapData> haps(total_hap_length[hap_buf]);
			add_haps(haps, batch.haps, i - cnt, cnt);
			haps_buffers.push_back(haps);

			haps_counts.push_back(cnt);
			cnt = 0;

			hap_buf++;

			total_hap_length.push_back(0);
		}

		cnt++;
		total_hap_length[hap_buf] += hap_length_padded;
	}

	inaccel::vector<HapData> haps(total_hap_length[hap_buf]);
	add_haps(haps, batch.haps, batch.num_haps - cnt, cnt);
	haps_buffers.push_back(haps);

	haps_counts.push_back(cnt);

	std::vector<inaccel::vector<ResultData>> results_buffers(reads_buffers.size() * haps_buffers.size());

	std::vector<::inaccel::session> sessions;

	for (int r = 0; r < reads_buffers.size(); r++) {
		for (int h = 0; h < haps_buffers.size(); h++) {
			int read_length = total_read_length[r];
			int hap_length = total_hap_length[h];
			int num_rows = max_read_length[r] + hap_length - 1;
			int num_results = reads_counts[r] * haps_counts[h];

			results_buffers[r * hap_buf + h] = inaccel::vector<ResultData>(reads_counts[r] * haps_counts[h]);

			// Create request for “pairhmm_driver” and “pairhmm_output” kernels
			inaccel::request pairhmm("com.wasai.genomics.pairhmm");

			pairhmm.arg(reads_buffers[r]).arg(haps_buffers[h]).arg(read_length).arg(hap_length).arg(num_rows).arg(results_buffers[r * hap_buf + h]).arg(num_results);

			// Submit the request
			sessions.push_back(inaccel::submit(pairhmm));

		}
	}

	int m = 0;
	// Wait for all requests to finish
	for(int req = 0; req < reads_buffers.size() * haps_buffers.size(); req++) {
		inaccel::wait(sessions[req]);
		// Save results
		for(int i = 0; i < results_buffers[req].size(); i++) {
			batch.results[m++] = results_buffers[req][i].result;
		}

	}


	return true;
}
