#include "Timer.hpp"

namespace MyFem {
  
#if(TIMERS_ON)
  
void TimerRegistry::addTimer(const Timer& timer) {
  std::lock_guard<std::mutex> lock(mutex_);
  auto it = timers_.find(timer.name_);
  if (it == timers_.end())
    timers_[timer.name_] = timer; // this does copy
  else {
    timers_[timer.name_].timerVal_secDouble_ += timer.timerVal_secDouble_; // direct use for efficiency?
  }
}

std::string TimerRegistry::timingReportStr() const {
  std::string out = "==== Timing Report ====\n";
  out += "Registry name: "+registryName_+"\n";
  double totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-totalStart_).count();
  out += "Total registry time | " + std::to_string(totalTime) + "s\n";
  for(const auto& [name, timer] : timers_) { // TODO: Can I somehow decide upon an order? Or just in general make it not unordered (but still efficient somehow?); Need to make a smart algorithm or something I think...
    double timerVal = timer.getTimerVal_secDouble();
    double shareOfTotal = timerVal/totalTime;
    out += "Timer "+name+" | "+std::to_string(timerVal)+"s | "+std::to_string(shareOfTotal*100)+"% | (started at "+std::to_string(std::chrono::duration<double>(timer.start_-totalStart_).count())+"s)\n"; // TODO: for things like this, we should make a general function which makes a table or soemthow spaces the things automatically (similar to matrix print)
  }

  return levelizeString(out, 1);
}

#else

std::string TimerRegistry::timingReportStr() const noexcept {
  std::string out = "==== Timing Report : Timers are off (TIMERSON == false) ====\n";
  out += "Registry name: "+registryName_+"\n";
  double totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-totalStart_).count();
  out += "Total registry time | " + std::to_string(totalTime) + "s\n";
  return levelizeString(out, 1);
}

#endif

}
