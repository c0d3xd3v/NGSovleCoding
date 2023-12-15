import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15

RowLayout {
    Layout.fillWidth: true
    Layout.preferredHeight: 36
    Layout.alignment: Qt.AlignCenter
    GridLayout {
        columnSpacing: 3
        rowSpacing: 0
        columns: 5
        Repeater {
            model: 3
            Rectangle {
                width: 32
                height: 32
                radius: 90
                clip:true
                Button {
                    id: control
                    flat: false
                    anchors.fill: parent
                    text: index
                    onClicked: {
                        view.currentIndex = index
                    }
                    /*
                    background: Rectangle {
                            radius: 90
                            implicitWidth: 60
                            implicitHeight: 60
                            opacity: enabled ? 1 : 0.3
                            color: control.down ? "#d0d0d0" : "#e0e0e0"
                        }
                    */
                }
            }
        }
    }
}
