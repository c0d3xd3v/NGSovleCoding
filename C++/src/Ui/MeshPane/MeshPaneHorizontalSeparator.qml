import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

RowLayout {
    Layout.alignment: Qt.AlignHCenter
    Layout.fillWidth: true
    Layout.preferredHeight: 1.0
    Item {
        Layout.preferredHeight: 1
        Layout.fillWidth: true
    }
    Rectangle{
        Layout.preferredHeight: 2.0
        Layout.fillWidth: true
        radius: 20
        color: "gray"
    }
    Item {
        Layout.preferredHeight: 1
        Layout.fillWidth: true
    }
}
