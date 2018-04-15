#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <ctime>

#include <hiredis/hiredis.h>
#include <hiredis/async.h>

#include "Redis.hh"
#include "ChannelData.hh"
#include "HeaderData.hh"
#include "Noise.hh"
#include "Analysis.h"
#include "ChannelMap.hh"
#include "FFT.hh"
#include "RedisData.hh"

using namespace daqAnalysis;
using namespace std;

Redis::Redis(std::vector<unsigned> &stream_take, std::vector<unsigned> &stream_expire, int snapshot_time, int waveform_input_size, bool timing): 
  _snapshot_time(snapshot_time),
  _stream_take(stream_take),
  _stream_expire(stream_expire),
  _stream_last(stream_take.size(), 0),
  _stream_send(stream_take.size(), false),
  _last_snapshot(0),

  // allocate and zero-initalize all of the metrics
  _rms(stream_take.size(), RedisRMS()),
  _baseline(stream_take.size(), RedisBaseline()),
  _dnoise(stream_take.size(), RedisDNoise()),
  _pulse_height(stream_take.size(), RedisPulseHeight()),
  _occupancy(stream_take.size(), RedisOccupancy()),

  _fem_scaled_sum_rms(stream_take.size(), StreamDataMean(ChannelMap::n_fem_per_board* ChannelMap::n_boards, 1)),

  // and the header stuff
  _frame_no(stream_take.size(), StreamDataMax(ChannelMap::n_fem_per_board* ChannelMap::n_boards)),
  _trigframe_no(stream_take.size(), StreamDataMax(ChannelMap::n_fem_per_board* ChannelMap::n_boards)),
  _event_no(stream_take.size(), StreamDataMax(ChannelMap::n_fem_per_board* ChannelMap::n_boards)),

  _fft_manager((waveform_input_size > 0) ? waveform_input_size: 0),
  _do_timing(timing)
{
  //context = redisAsyncConnect("127.0.0.1", 6379);
  context = redisConnect("127.0.0.1", 6379);
  std::cout << "REDIS CONTEXT: " << context << std::endl;
  if (context != NULL && context->err) {
    std::cerr << "Redis error: " <<  context->errstr << std::endl;
    exit(1);
  }
  _start = 0;
}

Redis::~Redis() {
  //redisAsyncDisconnect(context);
  redisFree(context);
}

void Redis::StartSend() {
  _now = std::time(nullptr);
  // first time startup
  if (_start == 0) _start = _now;

  for (size_t i = 0; i < _stream_send.size(); i++) {
    // calculate whether each stream is sending to redis on this event
    _stream_send[i] = _stream_last[i] != _now && (_now - _start) % _stream_take[i] == 0;
    if (_stream_send[i]) {
      _stream_last[i] = _now;
    }
  }
}

void Redis::FinishSend() {
  if (_do_timing) {
    _timing.Print();
  }
}

void Redis::SendHeaderData(vector<HeaderData> *header_data) {
  for (auto &header: *header_data) {
    // update header info in each stream
    for (size_t i = 0; i < _stream_take.size(); i++) {
      // index into the fem data cache
      unsigned fem_ind = header.Ind();
      _event_no[i].Add(fem_ind, header.event_number);
      _frame_no[i].Add(fem_ind, header.frame_number);
      _trigframe_no[i].Add(fem_ind, header.trig_frame_number);
    }
  }
  for (size_t i = 0; i < _stream_take.size(); i++) {
    // send headers to redis if need be
    if (_stream_send[i]) {
      SendHeader(i);
    }
  }
}

void Redis::SendHeader(unsigned stream_index) {
  void *reply; 
  if (_do_timing) {
    _timing.StartTime();
  }
  for (size_t fem_ind = 0; fem_ind < ChannelMap::n_fem_per_board* ChannelMap::n_boards; fem_ind++) {
    // TODO @INSTALLATION: implement translation from fem_ind to fem/board
    // TEMPORARY IMPLEMENTATION
    unsigned fem = fem_ind % ChannelMap::n_fem_per_board;
    unsigned board = fem_ind / ChannelMap::n_fem_per_board;

    reply = redisCommand(context, "SET stream/%i:%i:frame_no:board:%i:fem:%i %f", 
      _stream_take[stream_index], _now/_stream_take[stream_index], fem, board, _frame_no[stream_index].Take(fem_ind));
    freeReplyObject(reply);
    
    reply = redisCommand(context, "EXPIRE stream/%i:%i:frame_no:board:%i:fem:%i %i", 
      _stream_take[stream_index], _now/_stream_take[stream_index], fem, board, _stream_expire[stream_index]);
    freeReplyObject(reply);
    

    reply = redisCommand(context, "SET stream/%i:%i:event_no:board:%i:fem:%i %f", 
      _stream_take[stream_index], _now/_stream_take[stream_index], fem, board, _event_no[stream_index].Take(fem_ind));
    freeReplyObject(reply);
    
    reply = redisCommand(context, "EXPIRE stream/%i:%i:event_no:board:%i:fem:%i %i", 
      _stream_take[stream_index], _now/_stream_take[stream_index], fem, board, _stream_expire[stream_index]);
    freeReplyObject(reply);
    

    reply = redisCommand(context, "SET stream/%i:%i:trigframe_no:board:%i:fem:%i %f", 
      _stream_take[stream_index], _now/_stream_take[stream_index], fem, board, _trigframe_no[stream_index].Take(fem_ind));
    freeReplyObject(reply);
    
    reply = redisCommand(context, "EXPIRE stream/%i:%i:trigframe_no:board:%i:fem:%i %i", 
      _stream_take[stream_index], _now/_stream_take[stream_index], fem, board, _stream_expire[stream_index]);
    freeReplyObject(reply);
    
    if (_do_timing) {
      _timing.EndTime(&_timing.send_header_data);
    }
  }
}

// TEMPORARY: Still needed to _ssum_rms
// may remove later?
void Redis::SendFem(unsigned stream_index) {
  // TODO: Report failures

  void *reply;
  if (_do_timing) {
    _timing.StartTime();
  }
  for (unsigned fem_ind = 0; fem_ind < _fem_scaled_sum_rms[stream_index].Size(); fem_ind++) {
    // TODO @INSTALLATION: implement translation from fem_ind to fem/board
    // TEMPORARY IMPLEMENTATION
    unsigned fem = fem_ind % ChannelMap::n_fem_per_board;
    unsigned board = fem_ind / ChannelMap::n_fem_per_board;
    reply = redisCommand(context, "SET stream/%i:%i:scaled_sum_rms:board:%i:fem:%i %f", 
      _stream_take[stream_index], _now/_stream_take[stream_index], board, fem, _fem_scaled_sum_rms[stream_index].Take(fem_ind));
    freeReplyObject(reply);
           
    reply = redisCommand(context, "EXPIRE stream/%i:%i:scaled_sum_rms:board:%i:fem:%i %i", 
      _stream_take[stream_index], _now/_stream_take[stream_index], board, fem, _stream_expire[stream_index]);
    freeReplyObject(reply);
  }
  if (_do_timing) {
    _timing.EndTime(&_timing.send_fem_data);
  }
}

inline void pushFFTDat(redisContext *context, unsigned channel, float re, float im) {
  float dat = re * im;
  redisAppendCommand(context, "RPUSH snapshot:fft:wire:%i %f", channel, dat);
}

void Redis::Snapshot(vector<ChannelData> *per_channel_data, vector<NoiseSample> *noise) {
  size_t n_commands = 0;

  // record the time for reference
  redisAppendCommand(context, "SET snapshot:time %i", _now);
  n_commands ++;

  // stuff per channel
  for (auto &channel: *per_channel_data) { 
    if (_do_timing) {
      _timing.StartTime();
    }
    // store the waveform and fft's
    // also delete old lists
    redisAppendCommand(context, "DEL snapshot:waveform:wire:%i", channel.channel_no);
    
    for (auto dat: channel.waveform) {
      redisAppendCommand(context, "RPUSH snapshot:waveform:wire:%i %i", channel.channel_no, dat);
    }
    n_commands += channel.waveform.size() + 1;

    if (_do_timing) {
      _timing.EndTime(&_timing.send_waveform);
    }

    if (_do_timing) {
      _timing.StartTime();
    }

    redisAppendCommand(context, "DEL snapshot:fft:wire:%i", channel.channel_no);
    // use already calculated FFT if there
    if (channel.fft_real.size() != 0) {
      for (size_t i = 0; i < channel.fft_real.size(); i++) {
        pushFFTDat(context, channel.channel_no, channel.fft_real[i], channel.fft_imag[i]);
      }
      n_commands += channel.fft_real.size() + 1;
    }
    // otherwise calculate it yourself
    else {
      // re-allocate if necessary
      if (_fft_manager.InputSize() != channel.waveform.size()) {
        _fft_manager.Set(channel.waveform.size());
      }
      // fill up the fft data in
      for (size_t i = 0; i < channel.waveform.size(); i++) {
        double *input = _fft_manager.InputAt(i);
        *input = (double) channel.waveform[i];
      }
      _fft_manager.Execute();
      int adc_fft_size = _fft_manager.OutputSize();
      for (int i = 0; i < adc_fft_size; i++) {
        pushFFTDat(context, channel.channel_no, _fft_manager.ReOutputAt(i), _fft_manager.ImOutputAt(i));
      }
      n_commands += adc_fft_size + 1;
    }
    if (_do_timing) {
      _timing.EndTime(&_timing.send_fft);
    }
  }

  if (_do_timing) {
    _timing.StartTime();
  }
  // also store the noise correlation matrix if taking a snapshot
  redisAppendCommand(context, "DEL snapshot:correlation");
  
  // Only calculate the upper-right half of the matrix since it is symmetric
  // The index 'k' into the list of the i-th sample with the j-th sample (where i <= j) is:
  // k = ((n+1)*n/2) - (n-i+1)*(n-i)/2 + j - i
  for (size_t i = 0; i < noise->size(); i++) {
    for (size_t j = i; j < noise->size(); j++) {
      float correlation = (*noise)[i].Correlation((*per_channel_data)[i].waveform, (*noise)[j], (*per_channel_data)[j].waveform);
      redisAppendCommand(context, "RPUSH snapshot:correlation %f",  correlation);
    }
  }
  n_commands += 1 +(noise->size() + 1) * noise->size() / 2;
  if (_do_timing) {
    _timing.EndTime(&_timing.correlation);
  }

  if (_do_timing) {
    _timing.StartTime();
  }
  FinishPipeline(n_commands);
  if (_do_timing) {
    _timing.EndTime(&_timing.clear_pipeline);
  }
  
}

void Redis::SendChannelData(vector<ChannelData> *per_channel_data, vector<NoiseSample> *noise_samples) {
  if (_do_timing) {
    _timing.StartTime();
  }

  unsigned n_channels = per_channel_data->size();
  unsigned n_fem = ChannelMap::NFEM();
  bool at_end_of_detector = false;

  // iterate over boards and fems
  for (unsigned board = 0; board < daqAnalysis::ChannelMap::n_boards; board++) {
    for (unsigned fem = 0; fem < daqAnalysis::ChannelMap::n_fem_per_board; fem++) {
      // TODO @INSTALLATION: Will this still work?
      // index into the fem data cache
      unsigned fem_ind = board * ChannelMap::n_fem_per_board + fem;

      // detect if at end of fem's
      if (fem_ind >= n_fem) {
        at_end_of_detector = true;
        break;
      }
      // keeping track of stuff to calculate sum rms
      std::vector<std::vector<int16_t> *> waveforms {};
      std::vector<daqAnalysis::NoiseSample *> this_noise_samples {};

      for (unsigned channel = 0; channel < daqAnalysis::ChannelMap::n_channel_per_fem; channel++) {
        // get the wire number
	daqAnalysis::ChannelMap::board_channel board_channel {board, fem, channel};
	uint16_t wire = daqAnalysis::ChannelMap::Channel2Wire(board_channel);
        // detect if at end of detector
        if (wire >= n_channels) {
          at_end_of_detector = true;
          break;
        }
 
        // fill metrics for each stream
        for (size_t i = 0; i < _stream_take.size(); i++) {
          _rms[i].Add((*per_channel_data)[wire], board, fem_ind, wire);
          _baseline[i].Add((*per_channel_data)[wire], board, fem_ind, wire);
          _dnoise[i].Add((*per_channel_data)[wire], board, fem_ind, wire);
          _pulse_height[i].Add((*per_channel_data)[wire], board, fem_ind, wire);
          _occupancy[i].Add((*per_channel_data)[wire], board, fem_ind, wire);
        }

        // collect waveforms for sum rms calculation
        waveforms.push_back(&(*per_channel_data)[wire].waveform);
        this_noise_samples.push_back(&(*noise_samples)[wire]);
      }

      // sum rms!
      float ssum_rms = NoiseSample::ScaledSumRMS(this_noise_samples, waveforms); 
      // and then clear out containers
      waveforms.clear();
      this_noise_samples.clear();


      // update all the streams for sum rms
      for (size_t i = 0; i < _stream_take.size(); i++) {
        _fem_scaled_sum_rms[i].Add(fem_ind, ssum_rms);
      }
      // if at end, break
      if (at_end_of_detector) {
        break;
      }
  
    }
    // if at end, break
    if (at_end_of_detector) {
      break;
    }
  }
  if (_do_timing) {
    _timing.EndTime(&_timing.copy_data);
  }

  for (size_t i = 0; i < _stream_take.size(); i++) {
    // Send stuff to redis if it's time
    if (_stream_send[i]) {
      // send the fem stuff (scaled sum rms)
      SendFem(i);
      if (_do_timing) {
        _timing.StartTime();
      }
      // metrics control the sending of everything else
      _rms[i].Send(context, _now, _stream_take[i], _stream_expire[i]);
      _baseline[i].Send(context, _now, _stream_take[i], _stream_expire[i]);
      _dnoise[i].Send(context, _now, _stream_take[i], _stream_expire[i]);
      _occupancy[i].Send(context, _now, _stream_take[i], _stream_expire[i]);
      _pulse_height[i].Send(context, _now, _stream_take[i], _stream_expire[i]);
      if (_do_timing) {
        _timing.EndTime(&_timing.send_metrics);
      }
    }

    // the metric was taken iff it was sent to redis
    bool taken = _stream_send[i];
    // update all of the metrics
    _rms[i].Update(taken);
    _baseline[i].Update(taken);
    _dnoise[i].Update(taken);
    _pulse_height[i].Update(taken);
    _occupancy[i].Update(taken);
    // and sum rms
    _fem_scaled_sum_rms[i].Update(taken);
  }

  if (_do_timing) {
    _timing.StartTime();
  }
  // snapshots
  bool take_snapshot = _snapshot_time > 0 && (_now - _start) % _snapshot_time == 0 && _last_snapshot != _now;
  if (take_snapshot) {
    Snapshot(per_channel_data, noise_samples);
    _last_snapshot = _now;
  }
}

// clear out a sequence of Append Commands
void Redis::FinishPipeline(size_t n_commands) {
  // TODO: Error Handling
  void *reply;
  for (size_t i =0; i <n_commands; i++) {
    redisGetReply(context, &reply); 
    freeReplyObject(reply);
  }
}

void RedisTiming::StartTime() {
  start = clock();
}
void RedisTiming::EndTime(float *field) {
  *field += float(clock() - start)/CLOCKS_PER_SEC;
}
void RedisTiming::Print() {
  std::cout << "COPY DATA: " << copy_data << std::endl;
  std::cout << "SEND METRICS: " << send_metrics << std::endl;
  std::cout << "SEND FEM    : " << send_fem_data << std::endl;
  std::cout << "SEND HEADER : " << send_header_data << std::endl;
  std::cout << "SEND WAVEFORM " << send_waveform << std::endl;
  std::cout << "SEND FFT    : " << send_fft << std::endl;
  std::cout << "CORRELATION : " << correlation << std::endl;
  std::cout << "CLEAR PIPE  : " << clear_pipeline << std::endl;
}
