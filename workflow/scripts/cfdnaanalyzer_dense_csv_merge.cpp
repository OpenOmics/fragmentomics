#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

class CsvReader {
 public:
  explicit CsvReader(const std::string& path) : path_(path), buffer_(1 << 20) {
    input_.rdbuf()->pubsetbuf(buffer_.data(), buffer_.size());
    input_.open(path, std::ios::binary);
    if (!input_) throw std::runtime_error("cannot open " + path);
  }

  bool next(std::string& field, bool& record_end) {
    field.clear();
    record_end = false;
    int first = input_.get();
    if (first == EOF) return false;
    if (first == '"') {
      while (true) {
        int c = input_.get();
        if (c == EOF) throw std::runtime_error("unterminated quoted field in " + path_);
        if (c == '"') {
          int following = input_.peek();
          if (following == '"') {
            input_.get();
            field.push_back('"');
            continue;
          }
          c = input_.get();
          if (c == ',') return true;
          if (c == '\n' || c == EOF) {
            record_end = true;
            return true;
          }
          if (c == '\r') {
            if (input_.peek() == '\n') input_.get();
            record_end = true;
            return true;
          }
          throw std::runtime_error("characters after quoted field in " + path_);
        }
        field.push_back(static_cast<char>(c));
      }
    }
    int c = first;
    while (true) {
      if (c == ',') return true;
      if (c == '\n') {
        record_end = true;
        return true;
      }
      if (c == '\r') {
        if (input_.peek() == '\n') input_.get();
        record_end = true;
        return true;
      }
      if (c == '"') throw std::runtime_error("quote in unquoted field in " + path_);
      field.push_back(static_cast<char>(c));
      c = input_.get();
      if (c == EOF) {
        record_end = true;
        return true;
      }
    }
  }

 private:
  std::string path_;
  std::vector<char> buffer_;
  std::ifstream input_;
};

void write_csv_field(std::ostream& output, const std::string& value) {
  if (value.find_first_of(",\"\r\n") == std::string::npos) {
    output << value;
    return;
  }
  output.put('"');
  for (char c : value) {
    if (c == '"') output.put('"');
    output.put(c);
  }
  output.put('"');
}

class HeaderIterator {
 public:
  explicit HeaderIterator(const std::string& path) : path_(path), reader_(path) {
    bool end = false;
    std::string field;
    if (!reader_.next(field, end) || end || field != "sample")
      throw std::runtime_error("header does not begin with sample in " + path);
    if (!reader_.next(field, end) || field != "label")
      throw std::runtime_error("second header field is not label in " + path);
    header_finished_ = end;
  }

  bool advance() {
    if (header_finished_) {
      valid_ = false;
      return false;
    }
    bool end = false;
    std::string raw;
    if (!reader_.next(raw, end)) throw std::runtime_error("unexpected EOF in header " + path_);
    if (raw.find_first_of("\r\n") != std::string::npos)
      throw std::runtime_error("newline in feature name in " + path_);
    if (have_raw_ && raw < previous_raw_)
      throw std::runtime_error("feature header is not sorted in " + path_ +
                               ": " + raw + " follows " + previous_raw_);
    if (have_raw_ && raw == previous_raw_) {
      ++duplicate_number_;
      current_ = raw + "." + std::to_string(duplicate_number_);
    } else {
      duplicate_number_ = 0;
      current_ = raw;
    }
    // Resolve collisions such as raw columns a,a,a.1 without a global set.
    if (have_emitted_ && current_ <= previous_emitted_) current_ += ".1";
    if (have_emitted_ && current_ <= previous_emitted_)
      throw std::runtime_error("duplicate-name normalization is not monotonic in " + path_);
    previous_raw_ = std::move(raw);
    previous_emitted_ = current_;
    have_raw_ = have_emitted_ = valid_ = true;
    header_finished_ = end;
    return true;
  }

  const std::string& value() const { return current_; }
  bool valid() const { return valid_; }
  CsvReader& reader() { return reader_; }

 private:
  std::string path_;
  CsvReader reader_;
  std::string current_, previous_raw_, previous_emitted_;
  std::size_t duplicate_number_ = 0;
  bool have_raw_ = false, have_emitted_ = false, header_finished_ = false, valid_ = false;
};

struct HeapItem {
  std::string value;
  std::size_t source;
};
struct HeapGreater {
  bool operator()(const HeapItem& left, const HeapItem& right) const {
    if (left.value != right.value) return left.value > right.value;
    return left.source > right.source;
  }
};

std::vector<std::pair<std::string, std::string>> read_list(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open source list " + path);
  std::vector<std::pair<std::string, std::string>> entries;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto tab = line.find('\t');
    if (tab == std::string::npos || tab == 0 || tab + 1 == line.size())
      throw std::runtime_error("invalid source-list row: " + line);
    entries.emplace_back(line.substr(0, tab), line.substr(tab + 1));
  }
  if (entries.empty()) throw std::runtime_error("source list is empty");
  return entries;
}

void schema_mode(const std::string& list_path, const std::string& schema_path,
                 const std::string& header_path) {
  auto entries = read_list(list_path);
  std::vector<std::unique_ptr<HeaderIterator>> headers;
  headers.reserve(entries.size());
  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapGreater> heap;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    headers.emplace_back(std::make_unique<HeaderIterator>(entries[i].second));
    if (headers.back()->advance()) heap.push({headers.back()->value(), i});
  }
  std::ofstream schema(schema_path, std::ios::binary);
  std::ofstream header(header_path, std::ios::binary);
  if (!schema || !header) throw std::runtime_error("cannot create schema outputs");
  header << "sample,label";
  std::size_t features = 0;
  std::string previous;
  bool have_previous = false;
  while (!heap.empty()) {
    HeapItem item = heap.top();
    heap.pop();
    if (!have_previous || item.value != previous) {
      schema << item.value << '\n';
      header.put(',');
      write_csv_field(header, item.value);
      previous = item.value;
      have_previous = true;
      ++features;
    }
    auto& source = *headers[item.source];
    if (source.advance()) heap.push({source.value(), item.source});
  }
  header.put('\n');
  schema.flush();
  header.flush();
  if (!schema || !header) throw std::runtime_error("failed while writing schema outputs");
  std::cout << "features=" << features << "\ncolumns=" << features + 2 << '\n';
}

void row_mode(const std::string& source_path, const std::string& expected_sample,
              const std::string& schema_path, const std::string& output_path) {
  HeaderIterator header(source_path);
  bool have_header = header.advance();
  std::ifstream schema(schema_path);
  if (!schema) throw std::runtime_error("cannot open schema " + schema_path);
  std::vector<std::uint8_t> present;
  present.reserve(16000000);
  std::string feature;
  while (std::getline(schema, feature)) {
    if (!feature.empty() && feature.back() == '\r') feature.pop_back();
    while (have_header && header.value() < feature)
      throw std::runtime_error("source feature is absent from schema: " + header.value());
    if (have_header && header.value() == feature) {
      present.push_back(1);
      have_header = header.advance();
    } else {
      present.push_back(0);
    }
  }
  if (have_header) throw std::runtime_error("schema ended before source header");

  bool ended = false;
  std::string sample, label, value;
  CsvReader& row = header.reader();
  if (!row.next(sample, ended)) {
    std::ofstream output(output_path, std::ios::binary);
    if (!output) throw std::runtime_error("cannot create row output " + output_path);
    std::cout << "features=" << present.size() << "\ncolumns=" << present.size() + 2
              << "\nsource_features=0\nrows=0\n";
    return;
  }
  if (ended)
    throw std::runtime_error("missing data row in " + source_path);
  if (sample != expected_sample)
    throw std::runtime_error("data-row sample differs in " + source_path);
  if (!row.next(label, ended)) throw std::runtime_error("missing label in " + source_path);

  std::ofstream output(output_path, std::ios::binary);
  if (!output) throw std::runtime_error("cannot create row output " + output_path);
  write_csv_field(output, sample);
  output.put(',');
  write_csv_field(output, label);
  std::size_t source_features = 0;
  for (std::uint8_t bit : present) {
    output.put(',');
    if (bit) {
      if (ended || !row.next(value, ended))
        throw std::runtime_error("data row ended before its header in " + source_path);
      write_csv_field(output, value);
      ++source_features;
    }
  }
  if (!ended) {
    if (row.next(value, ended))
      throw std::runtime_error("data row has more fields than its header in " + source_path);
  }
  output.put('\n');
  output.flush();
  if (!output) throw std::runtime_error("failed while writing row output");
  std::cout << "features=" << present.size() << "\ncolumns=" << present.size() + 2
            << "\nsource_features=" << source_features << "\nrows=1\n";
}

std::string argument(int argc, char** argv, const std::string& name) {
  for (int i = 2; i + 1 < argc; ++i)
    if (argv[i] == name) return argv[i + 1];
  throw std::runtime_error("missing argument " + name);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) throw std::runtime_error("usage: helper schema|row [arguments]");
    const std::string mode = argv[1];
    if (mode == "schema") {
      schema_mode(argument(argc, argv, "--list"), argument(argc, argv, "--schema"),
                  argument(argc, argv, "--header"));
    } else if (mode == "row") {
      row_mode(argument(argc, argv, "--source"), argument(argc, argv, "--sample"),
               argument(argc, argv, "--schema"), argument(argc, argv, "--output"));
    } else {
      throw std::runtime_error("unknown mode: " + mode);
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
